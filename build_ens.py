#!/usr/bin/env python3
"""
build_ens.py - 把 5 個集成預報的熱帶氣旋路徑轉成一個細 JSON (ens.json)

來源 (全部公開、免登入):
  ECMWF IFS ENS    https://data.ecmwf.int/forecasts   (BUFR, CC-BY-4.0, ECMWF)
  ECMWF AIFS ENS   https://data.ecmwf.int/forecasts   (BUFR, CC-BY-4.0, ECMWF)
  NOAA GEFS        https://nomads.ncep.noaa.gov/pub/data/nccf/com/ens_tracker/prod/gefs.*
  NOAA AIGEFS      https://nomads.ncep.noaa.gov/pub/data/nccf/com/ens_tracker/prod/aigefs.*
  Google FNV3      https://deepmind.google.com/science/weatherlab/download/cyclones/FNV3/ensemble/...
                   (Google 實驗數據,有使用條款,見 README)

做法: 以 ECMWF 嘅有編號風暴 (編號 < 70, 例如 33W/34W/35W) 做風暴名單,
其他來源用「初始位置相差 < 3 度」去配對 (唔用盆地碼, 因為 NOLO 喺 NOAA 係 EP15)。
配唔到嘅風暴,嗰個模型就唔畫。冇任何估算或補數據。
"""
import json, math, sys, os, re, datetime as dt, urllib.request, urllib.error, tempfile

UA = {"User-Agent": "Mozilla/5.0 (ens-track-builder)"}
STEP_H = 6            # 輸出每 6 小時一點
BOX = (0, 50, 95, 200) # lat min,max, lon min,max (東經, 西經已轉成 +360)

def http(url, binary=False, timeout=60):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        d = r.read()
    return d if binary else d.decode("utf-8", "replace")

def head_ok(url):
    try:
        req = urllib.request.Request(url, headers=UA, method="HEAD")
        urllib.request.urlopen(req, timeout=30).close()
        return True
    except Exception:
        return False

def norm_lon(lon):
    return lon + 360 if lon < -120 else lon

def runs_newest_first(now=None, days=2):
    now = now or dt.datetime.utcnow()
    out = []
    for d in range(days + 1):
        day = (now - dt.timedelta(days=d)).date()
        for hh in (18, 12, 6, 0):
            t = dt.datetime(day.year, day.month, day.day, hh)
            if t <= now:
                out.append(t)
    return out

# ---------------- ECMWF (BUFR) ----------------
def ecmwf_fetch(stream):  # stream = "ifs" or "aifs-ens"
    for t in runs_newest_first():
        ymd, hh = t.strftime("%Y%m%d"), t.strftime("%H")
        for step in ("360h", "144h"):
            url = f"https://data.ecmwf.int/forecasts/{ymd}/{hh}z/{stream}/0p25/enfo/{ymd}{hh}0000-{step}-enfo-tf.bufr"
            if head_ok(url):
                return t, url, http(url, binary=True, timeout=120)
    raise RuntimeError("ECMWF " + stream + " 搵唔到最近嘅 tf.bufr")

def ecmwf_parse(blob):
    import eccodes as ec
    path = tempfile.mktemp(suffix=".bufr")
    open(path, "wb").write(blob)
    storms = {}
    with open(path, "rb") as fh:
        while True:
            h = ec.codes_bufr_new_from_file(fh)
            if h is None:
                break
            ec.codes_set(h, "unpack", 1)
            sid = str(ec.codes_get_array(h, "stormIdentifier")[0]).strip()
            m = re.match(r"^(\d+)([A-Z])$", sid)
            if not m or int(m.group(1)) >= 70:
                ec.codes_release(h); continue
            name = str(ec.codes_get_array(h, "longStormName")[0]).strip()
            ns = ec.codes_get(h, "numberOfSubsets")
            def arr(key):
                try:
                    return ec.codes_get_array(h, key)
                except Exception:
                    return None
            K = 0
            ent = []   # 每個 #k# 位置: (lat陣列, lon陣列, 重要性陣列)
            while True:
                la = arr(f"#{K+1}#latitude")
                if la is None:
                    break
                K += 1
                ent.append((la, arr(f"#{K}#longitude"), arr(f"#{K}#meteorologicalAttributeSignificance")))
            tp = []
            j = 1
            while True:
                a_ = arr(f"#{j}#timePeriod")
                if a_ is None:
                    break
                tp.append(int(a_[0])); j += 1
            hours = [0] + tp
            def pick(a_, i):
                return a_[i] if len(a_) > 1 else a_[0]
            members = []
            for s in range(ns):
                pts = []
                ci = 0
                for la, lo, sg in ent:
                    if int(pick(sg, s)) != 1:
                        continue
                    if ci < len(hours):
                        hr = hours[ci]
                        v1, v2 = float(pick(la, s)), float(pick(lo, s))
                        if abs(v1) < 1e5 and abs(v2) < 1e5 and hr % STEP_H == 0:
                            pts.append([hr, round(v1, 1), round(norm_lon(v2), 1)])
                    ci += 1
                if pts:
                    members.append(pts)
            y = ec.codes_get_array(h, "year")[0]; mo = ec.codes_get_array(h, "month")[0]
            d = ec.codes_get_array(h, "day")[0]; hr = ec.codes_get_array(h, "hour")[0]
            storms[sid] = {"name": name, "members": members, "init": [members[0][0][1], members[0][0][2]] if members else None}
            ec.codes_release(h)
    os.remove(path)
    return storms

# ---------------- ATCF (NOAA / Google) ----------------
def lat_s(s): v = int(s[:-1]) / 10; return v if s[-1] == "N" else -v
def lon_s(s): v = int(s[:-1]) / 10; return v if s[-1] == "E" else -v

def atcf_parse(text, run_ymdh):
    """回傳 {(basin,num): {member_tag: {tau:(lat,lon)}}}"""
    out = {}
    for line in text.splitlines():
        if line.startswith("#"):
            continue
        p = [x.strip() for x in line.split(",")]
        if len(p) < 9 or p[2] != run_ymdh:
            continue
        try:
            tau = int(p[5]); la = lat_s(p[6]); lo = norm_lon(lon_s(p[7]))
        except Exception:
            continue
        out.setdefault((p[0], p[1]), {}).setdefault(p[4], {})[tau] = (la, lo)
    return out

def atcf_to_storms(parsed):
    res = {}
    for key, mem in parsed.items():
        members = []
        for tag, taus in sorted(mem.items()):
            if 0 not in taus:
                continue
            pts = [[t, round(v[0], 1), round(v[1], 1)] for t, v in sorted(taus.items()) if t % STEP_H == 0]
            if pts:
                members.append(pts)
        if members:
            res[key] = {"members": members, "init": [members[0][0][1], members[0][0][2]]}
    return res

def noaa_fetch(kind, prefix_files):
    base = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/ens_tracker/prod"
    for t in runs_newest_first():
        ymd, hh = t.strftime("%Y%m%d"), t.strftime("%H")
        probe = f"{base}/{kind}.{ymd}/{hh}/tctrack/{prefix_files[0]}.t{hh}z.cyclone.trackatcfunix"
        if not head_ok(probe):
            continue
        texts = []
        for pf in prefix_files:
            try:
                texts.append(http(f"{base}/{kind}.{ymd}/{hh}/tctrack/{pf}.t{hh}z.cyclone.trackatcfunix"))
            except Exception:
                pass
        return t, "\n".join(texts)
    raise RuntimeError("NOAA " + kind + " 搵唔到最近 run")

def google_fetch():
    base = "https://deepmind.google.com/science/weatherlab/download/cyclones/FNV3/ensemble/paired/atcf"
    for t in runs_newest_first():
        if t.hour not in (0, 6, 12, 18):
            continue
        url = f"{base}/FNV3_{t:%Y_%m_%dT%H}_00_atcf_a_deck.txt"
        try:
            return t, http(url, timeout=120)
        except urllib.error.HTTPError:
            continue
    raise RuntimeError("Google FNV3 搵唔到最近 run")

# ---------------- 主程式 ----------------
def dist(a, b):
    return math.hypot(a[0] - b[0], (a[1] - b[1]) * math.cos(math.radians(a[0])))

def match(master, storms, maxdeg=3.0):
    """master: {id:{name,init}} ; storms: {key:{members,init}} -> {master_id: storm}"""
    out = {}
    for mid, ms in master.items():
        best, bd = None, 1e9
        for k, st in storms.items():
            d = dist(ms["init"], st["init"])
            if d < bd:
                best, bd = k, d
        if best is not None and bd <= maxdeg:
            out[mid] = storms[best]
    return out

def main():
    models = []
    only = set(sys.argv[1:])  # 可指定只做某幾個: ecmwf aifs gefs aigefs fnv3
    def want(k): return not only or k in only

    master = None
    # ECMWF IFS ENS
    ecm = {}
    if want("ecmwf") or want("aifs") or not only:
        pass
    results = {}
    if want("ecmwf"):
        t, url, blob = ecmwf_fetch("ifs")
        st = ecmwf_parse(blob); results["ecmwf"] = (t, st, url)
    if want("aifs"):
        t, url, blob = ecmwf_fetch("aifs-ens")
        st = ecmwf_parse(blob); results["aifs"] = (t, st, url)

    # 風暴名單: 以 ECMWF (優先 IFS ENS) 做標準
    ref = (results.get("ecmwf") or results.get("aifs"))
    if not ref:
        raise SystemExit("需要 ecmwf 或 aifs 做風暴名單")
    master = {k: v for k, v in ref[1].items()
              if v["init"] and BOX[0] <= v["init"][0] <= BOX[1] and BOX[2] <= v["init"][1] <= BOX[3]}

    def add(mid, label, color, run, storms, src, lic):
        ms = match(master, storms) if mid not in ("ecmwf", "aifs") else {k: storms[k] for k in master if k in storms}
        models.append({
            "id": mid, "label": label, "color": color, "run": run.strftime("%Y-%m-%dT%H:00Z"),
            "source": src, "license": lic,
            "storms": {master[k]["name"]: {"members": v["members"]} for k, v in ms.items()},
        })

    if "ecmwf" in results:
        t, st, url = results["ecmwf"]
        add("ecmwf", "ECMWF ENS", "#38bdf8", t, st, "ECMWF open data", "CC-BY-4.0 ECMWF")
    if "aifs" in results:
        t, st, url = results["aifs"]
        add("aifs", "ECMWF AIFS ENS", "#a78bfa", t, st, "ECMWF open data", "CC-BY-4.0 ECMWF")
    if want("gefs"):
        t, txt = noaa_fetch("gefs", [f"ap{n:02d}" for n in range(1, 31)] + ["ac00"])
        add("gefs", "GEFS", "#34d399", t, atcf_to_storms(atcf_parse(txt, t.strftime("%Y%m%d%H"))), "NOAA NCEP ens_tracker", "NOAA public data")
    if want("aigefs"):
        t, txt = noaa_fetch("aigefs", [f"a{n:03d}" for n in range(0, 31)])
        add("aigefs", "AIGEFS", "#fbbf24", t, atcf_to_storms(atcf_parse(txt, t.strftime("%Y%m%d%H"))), "NOAA NCEP ens_tracker", "NOAA public data")
    if want("fnv3"):
        t, txt = google_fetch()
        add("fnv3", "WN3C (FNV3)", "#f472b6", t, atcf_to_storms(atcf_parse(txt, t.strftime("%Y%m%d%H"))), "Google DeepMind Weather Lab", "GDM Real-Time Weather Forecasting Experimental Data Terms of Use")

    out = {
        "generated": dt.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
        "stormNames": sorted({v["name"] for v in master.values()}),
        "models": models,
    }
    path = os.environ.get("ENS_OUT", "ens.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, separators=(",", ":"))
    print("wrote", path, os.path.getsize(path), "bytes")
    for m in models:
        print(m["id"], m["run"], {n: len(s["members"]) for n, s in m["storms"].items()})

if __name__ == "__main__":
    main()
