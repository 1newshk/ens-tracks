import unittest,datetime as dt
from build_genesis import cluster,distance
R=dt.datetime(2026,10,10,6,tzinfo=dt.timezone.utc)
def event(m,t,lat=15,lon=140,id='70W'):
 p={'tau':t,'lat':lat,'lon':lon,'validTime':(R+dt.timedelta(hours=t)).isoformat()}
 return {'member':m,'featureId':f'{m}:{id}:{t}','sourceStormId':id,'firstDetection':p,'track':[p]}
class Tests(unittest.TestCase):
 def test_same_local_id_not_identity(self):
  a=cluster([event(1,24),event(2,24,35)],{'complete':True},R);self.assertEqual(len(a),2)
 def test_same_member_never_double_counted(self):
  a=cluster([event(1,24),event(1,30,id='71W'),event(2,24)],{'complete':True},R);self.assertTrue(all(len({f['member'] for f in c['events']})==c['memberCount'] for c in a))
 def test_windows_bins(self):
  a=cluster([event(1,48),event(2,54),event(3,90)],{'complete':True},R)[0];self.assertEqual(a['counts'],{'0_48':1,'0_120':3});self.assertEqual(a['timeBins'],[1,1,1]);self.assertEqual(a['probability']['0_120'],.06)
 def test_missing_not_zero(self):
  a=cluster([event(1,24)],{'complete':False},R)[0];self.assertIsNone(a['probability']['0_48']);self.assertIsNone(a['notDetectedWithin120'])
 def test_no_transitive_chains(self):
  a=cluster([event(1,24,15,140),event(2,24,15,145),event(3,24,15,150)],{'complete':True},R);self.assertEqual(len(a),2)
 def test_dateline(self):
  self.assertLess(distance({'lat':0,'lon':179},{'lat':0,'lon':-179}),230)
 def test_time_split(self):self.assertEqual(len(cluster([event(1,12),event(2,120)],{'complete':True},R)),2)
if __name__=='__main__':unittest.main()
