"""Compare retained same-mesh outputs without reconstructing any equations."""
import csv,json,math
from study import BASE,dump
OUT=BASE/'gpu-v2'
def compare(a,b):
 def read(p):
  with p.open(encoding='utf-8-sig') as f:return list(csv.DictReader(f))
 x=read(a);y=read(b);assert len(x)==len(y)
 assert all(i['cell']==j['cell'] and i['area']==j['area'] and i['x']==j['x'] and i['y']==j['y'] for i,j in zip(x,y))
 dp=[float(j['p'])-float(i['p']) for i,j in zip(x,y)]
 area=sum(float(i['area']) for i in x);gauge=sum(float(i['area'])*d for i,d in zip(x,dp))/area
 return dict(cells=len(x),maxVelocityDifferenceOverUref=max(math.hypot(float(i['u'])-float(j['u']),float(i['v'])-float(j['v']))/.5 for i,j in zip(x,y)),pressureConstantDifference=gauge,maxPressureDifferenceAfterConstant=max(abs(d-gauge) for d in dp),pressureUnits='m2/s2',sources=[str(a),str(b)])
def main():
 pairs=[('same65k5steps','gmsh65k-simple-ic0-cpu-n5-min10','gmsh65k-simple-ic0-torch-cuda-n5-min10'),('full9kCoupledCpuGpu','shape160-coupled-ic0-cpu-n100-min1','shape160-coupled-ic0-torch-cuda-n100-min1'),('minimumIterations','shape160-coupled-aggregation-cpu-n100-min10','shape160-coupled-aggregation-cpu-n100-min1')]
 result={}
 for key,a,b in pairs:
  x=OUT/(a+'.cells.csv');y=OUT/(b+'.cells.csv')
  if x.exists() and y.exists():result[key]=compare(x,y)
 a=BASE/'solutions/graded128-simple-strict-ic0-1e-6.cells.csv';b=OUT/'graded128-simple-aggregation-cpu-n6000-min10.cells.csv'
 if b.exists():result['simplePreconditioner']=compare(a,b)
 dump(OUT/'field-differences.json',result);print(json.dumps(result,indent=2))
if __name__=='__main__':main()
