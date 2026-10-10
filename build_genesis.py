#!/usr/bin/env python3
"""Separate ECMWF genesis output. Does not alter build_ens.py or ens.json.
First identification is a tracker event, NOT official formation or a calibrated probability.
Supported BUFR layout verified against IFS/AIFS 2026-10-10 06Z; fails closed on changes.
"""
import argparse,datetime as dt,json,math,os,re,tempfile,urllib.request,hashlib
from pathlib import Path
UTC=dt.timezone.utc
EXPECTED=set(range(1,51));BOX=(0,50,95,200)
def distance(a,b):
 p,q=map(math.radians,[a['lat'],b['lat']]);dl=math.radians((a['lon']-b['lon']+180)%360-180)
 return 6371*2*math.asin(min(1,math.sqrt(math.sin((p-q)/2)**2+math.cos(p)*math.cos(q)*math.sin(dl/2)**2)))
def valid(x):return math.isfinite(float(x)) and abs(float(x))<1e5
def fetch(model):
 now=dt.datetime.now(UTC)
 for back in range(13):
  t=now.replace(hour=(now.hour//6)*6,minute=0,second=0,microsecond=0)-dt.timedelta(hours=back*6)
  for horizon in (360,144):
   url=f'https://data.ecmwf.int/forecasts/{t:%Y%m%d}/{t:%H}z/{model}/0p25/enfo/{t:%Y%m%d%H}0000-{horizon}h-enfo-tf.bufr'
   try:
    with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'genesis-builder/1.0'}),timeout=100) as r:
     blob=r.read();length=r.headers.get('Content-Length')
    if length and len(blob)!=int(length):raise ValueError('partial download')
    if not blob.startswith(b'BUFR') or not blob.endswith(b'7777'):raise ValueError('invalid BUFR envelope')
    return t,url,blob,horizon
   except Exception:continue
 raise RuntimeError('no valid recent ECMWF TC BUFR')
def parse(blob,run):
 import eccodes as ec
 features=[];inventory=set();types={};message_count=0;errors=[];horizon=0
 with tempfile.NamedTemporaryFile() as tmp:
  tmp.write(blob);tmp.flush()
  with open(tmp.name,'rb') as f:
   while True:
    h=ec.codes_bufr_new_from_file(f)
    if h is None:break
    message_count+=1
    try:
     ec.codes_set(h,'unpack',1)
     def arr(k):return list(ec.codes_get_array(h,k))
     def pick(a,i):
      if len(a) not in (1,ns):raise ValueError('unexpected subset shape '+str(len(a)))
      return a[0] if len(a)==1 else a[i]
     ns=int(ec.codes_get(h,'numberOfSubsets'));members=arr('ensembleMemberNumber');ft=arr('ensembleForecastType')
     if len(members)==1:members=members*ns
     if len(ft)==1:ft=ft*ns
     if len(members)!=ns or len(ft)!=ns:raise ValueError('missing per-subset member metadata')
     year,month,day,hour=[int(arr(k)[0]) for k in ('year','month','day','hour')]
     if dt.datetime(year,month,day,hour,tzinfo=UTC)!=run:raise ValueError('mixed run')
     for m,t in zip(members,ft):
      if int(t)==4:inventory.add(int(m))
      types[int(m)]=int(t)
     sid=str(arr('stormIdentifier')[0]).strip();match=re.fullmatch(r'(\d+)([A-Z])',sid)
     if not match or not 70<=int(match[1])<=99:continue
     # The supported template has observed centre rank1, model-analysis centre rank2,
     # wind-max rank3, then forecast-centre/wind-max pairs. Do NOT use flattened latitude.
     analysis_lat=arr('#2#latitude');analysis_lon=arr('#2#longitude');analysis_sg=arr('#2#meteorologicalAttributeSignificance')
     times=[];k=1
     while True:
      try:time=arr(f'#{k}#timePeriod')
      except ec.CodesInternalError:break
      times.append(int(time[0]));k+=1
     if any(t<=0 for t in times) or times!=sorted(set(times)):raise ValueError('unsupported forecast times')
     horizon=max(horizon,max(times,default=0));centres=[]
     for j,tau in enumerate(times):
      rank=4+j*2;la=arr(f'#{rank}#latitude');lo=arr(f'#{rank}#longitude');sg=arr(f'#{rank}#meteorologicalAttributeSignificance')
      if any(int(v) not in (1,ec.CODES_MISSING_LONG) for v in sg):raise ValueError('unsupported centre significance')
      centres.append((tau,la,lo,sg))
     for i,(member,typ) in enumerate(zip(members,ft)):
      member=int(member)
      if int(typ)!=4 or member not in EXPECTED:continue
      # Exclude model features already identified at initialization, whether observed or not.
      init_la,init_lo=pick(analysis_lat,i),pick(analysis_lon,i)
      if valid(init_la) and valid(init_lo):continue
      track=[]
      for tau,la,lo,sg in centres:
       y,x=pick(la,i),pick(lo,i)
       if int(pick(sg,i))==1 and valid(y) and valid(x) and -90<=float(y)<=90 and -180<=float(x)<=360:
        track.append({'tau':tau,'validTime':(run+dt.timedelta(hours=tau)).isoformat(),'lat':float(y),'lon':float(x)%360})
      if not track:continue
      first=track[0]
      if not (BOX[0]<=first['lat']<=BOX[1] and BOX[2]<=first['lon']<=BOX[3]):continue
      features.append({'featureId':f'{run:%Y%m%d%H}:{member}:{sid}:{message_count}','sourceStormId':sid,'member':member,'memberType':'perturbed','category':'forecast_new','firstDetection':first,'track':track})
    except Exception as ex:errors.append({'message':message_count,'error':str(ex)[:150]})
    finally:ec.codes_release(h)
 complete=inventory==EXPECTED and not errors and message_count>0 and horizon>=120
 return features,{'expectedMembers':50,'availableMembers':len(inventory&EXPECTED),'memberIds':sorted(inventory&EXPECTED),'complete':complete,'status':'complete' if complete else 'incomplete','basis':'all 50 perturbed member subset IDs represented; successful full-file decode; >=120h run','decodeErrors':errors,'messages':message_count,'decodedMaxHour':horizon}
def cluster(features,coverage,run):
 # Deterministic complete-link greedy grouping avoids arbitrarily long neighbour chains.
 # Member-local IDs are never association keys. Nearby events may still split/merge; heuristic.
 groups=[]
 for f in sorted(features,key=lambda f:(f['firstDetection']['tau'],f['firstDetection']['lat'],f['firstDetection']['lon'],f['member'],f['featureId'])):
  candidates=[]
  for i,g in enumerate(groups):
   if any(x['member']==f['member'] for x in g):continue
   costs=[(distance(f['firstDetection'],x['firstDetection']),abs(f['firstDetection']['tau']-x['firstDetection']['tau'])) for x in g]
   if all(d<=600 and t<=48 for d,t in costs):candidates.append((sum(d for d,t in costs)/len(costs),i))
  if candidates:groups[min(candidates)[1]].append(f)
  else:groups.append([f])
 out=[]
 for i,g in enumerate(groups):
  # Per-member earliest event; safeguards counting even if grouping evolves.
  unique={}
  for f in g:
   if f['member'] not in unique or f['firstDetection']['tau']<unique[f['member']]['firstDetection']['tau']:unique[f['member']]=f
  g=list(unique.values());taus=[f['firstDetection']['tau'] for f in g]
  counts={'0_48':sum(t<=48 for t in taus),'0_120':sum(t<=120 for t in taus)}
  bins=[sum(0<t<=48 for t in taus),sum(48<t<=72 for t in taus),sum(72<t<=120 for t in taus)]
  out.append({'clusterId':f'{run:%Y%m%d%H}:C{i+1:03}','displayLabel':f'預測區{i+1}','quality':'heuristic_unvalidated','firstDetectionTimeRange':[min(taus),max(taus)],'spatialExtent':{'latMin':min(f['firstDetection']['lat'] for f in g),'latMax':max(f['firstDetection']['lat'] for f in g),'lonMin':min(f['firstDetection']['lon'] for f in g),'lonMax':max(f['firstDetection']['lon'] for f in g)},'centre':{'lat':sum(f['firstDetection']['lat'] for f in g)/len(g),'lon':sum(f['firstDetection']['lon'] for f in g)/len(g)},'events':g,'memberCount':len(g),'denominator':50,'counts':counts,'probability':{k:v/50 if coverage['complete'] else None for k,v in counts.items()},'timeBins':bins,'binEdgesHours':[0,48,72,120],'notDetectedWithin120':50-counts['0_120'] if coverage['complete'] else None})
 return sorted(out,key=lambda c:(-c['counts']['0_120'],c['clusterId']))
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--ifs-file');ap.add_argument('--aifs-file');ap.add_argument('--run');ap.add_argument('--out',default='genesis.json');a=ap.parse_args();models=[]
 for model,path in [('ifs',a.ifs_file),('aifs-ens',a.aifs_file)]:
  try:
   if path:
    run=dt.datetime.fromisoformat(a.run.replace('Z','+00:00'));blob=Path(path).read_bytes();url='local-validation-sample';requested=None
   else:run,url,blob,requested=fetch(model)
   features,cov=parse(blob,run);models.append({'id':model,'label':'歐洲物理集成' if model=='ifs' else '歐洲人工智能集成','run':run.isoformat(),'sourceUrl':url,'sha256':hashlib.sha256(blob).hexdigest(),'requestedHorizon':requested,'coverage':cov,'clusters':cluster(features,cov,run)})
  except Exception as ex:models.append({'id':model,'label':'歐洲物理集成' if model=='ifs' else '歐洲人工智能集成','run':None,'coverage':{'expectedMembers':50,'availableMembers':0,'complete':False,'status':'download_or_decode_failed','error':str(ex)[:200]},'clusters':[]})
 out={'schema':'genesis.v1','generated':dt.datetime.now(UTC).isoformat(),'simulation':False,'definition':'first valid future forecast centre for ECMWF genesis-category feature absent in model analysis; NOT official formation time','method':{'distanceKm':600,'timeToleranceHours':48,'algorithm':'deterministic complete-link greedy, at most one event per member per cluster','probability':'raw unique perturbed-member proportion, not calibrated; per model only'},'unsupportedModels':[{'label':'美國人工智能集成','status':'未提供／未核實'},{'label':'谷歌人工智能集成','status':'未提供／未核實'}],'models':models}
 path=Path(a.out);path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(out,ensure_ascii=False,separators=(',',':')),encoding='utf8');os.replace(tmp,path)
 print(json.dumps({'output':str(path),'models':[{'id':m['id'],'coverage':m['coverage'],'clusters':len(m['clusters'])} for m in models]},ensure_ascii=False))
if __name__=='__main__':main()
