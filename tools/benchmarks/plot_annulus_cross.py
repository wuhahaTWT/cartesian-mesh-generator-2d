"""Plot actual accepted annulus grids/fields and the declared analytic target."""
from pathlib import Path
import argparse
import sys,csv,json,math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection,LineCollection
from matplotlib.font_manager import FontProperties
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
root=Path(__file__).resolve().parents[2];study=root/'outputs/annulus-cross';out=args.output
out.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(root/'tools/benchmarks'));import annulus_cross as model
sys.path.insert(0,str(root/'tools/flow'));import native_mesh as nm
font=FontProperties(fname='/Library/Fonts/Arial Unicode.ttf')
plt.rcParams.update({'font.family':font.get_name(),'axes.unicode_minus':False,'font.size':12,'figure.facecolor':'white','axes.spines.top':False,'axes.spines.right':False})
r=np.linspace(.5,1,400);v=(1/r-r)/3;p=(r*r-1/r**2)/18-2/9*np.log(r)
fig,axes=plt.subplots(1,3,figsize=(15,4.9),layout='constrained')
x=np.linspace(-1.02,1.02,260);X,Y=np.meshgrid(x,x);R=np.hypot(X,Y);speed=(1/np.maximum(R,1e-9)-R)/3;speed=np.ma.masked_where((R<.5)|(R>1),speed)
image=axes[0].pcolormesh(X,Y,speed,cmap='viridis',vmin=0,vmax=.5,shading='auto',rasterized=True)
for rad in (.58,.72,.88):
 t=np.linspace(0,2*np.pi,400);axes[0].plot(rad*np.cos(t),rad*np.sin(t),color='white',alpha=.65,lw=.8)
 theta=math.pi/3;pos=(rad*math.cos(theta),rad*math.sin(theta));end=(rad*math.cos(theta+.16),rad*math.sin(theta+.16));axes[0].annotate('',xy=end,xytext=pos,arrowprops=dict(arrowstyle='-|>',color='white',lw=1.5))
axes[0].set_aspect('equal');axes[0].set_xlabel('x / m');axes[0].set_ylabel('y / m');axes[0].set_title('真圆解析流场：流线为同心圆',fontproperties=font)
fig.colorbar(image,ax=axes[0],label='速率 / (m/s)',fraction=.04,pad=.025)
axes[1].plot(r,v,color='#2474a8',lw=2.4);axes[1].scatter([.5,1],[.5,0],color='#2474a8',s=30)
axes[1].set_xlabel('半径 r / m');axes[1].set_ylabel('切向速度 / (m/s)');axes[1].set_title('内壁 0.5 m/s，外壁 0',fontproperties=font);axes[1].grid(alpha=.18)
axes[1].text(.69,.4,r'$u_\theta(r)=\frac{1}{3}(\frac{1}{r}-r)$',fontsize=18)
axes[2].plot(r,p,color='#b26a32',lw=2.4);axes[2].scatter([.5,1],[p[0],0],color='#b26a32',s=30)
axes[2].set_xlabel('半径 r / m');axes[2].set_ylabel('运动学压力 p/ρ / (m²/s²)');axes[2].set_title('压力向外升高；取外壁压力为 0',fontproperties=font);axes[2].grid(alpha=.18)
axes[2].annotate(f'内外压差：{abs(p[0]):.6f} m²/s²',xy=(.5,p[0]),xytext=(.58,-.044),fontproperties=font,fontsize=11,arrowprops=dict(arrowstyle='->',color='#b26a32'))
fig.suptitle('理想二维圆形 Couette 流：Ri = 0.5 m，Ro = 1 m，ν = 0.1 m²/s',fontproperties=font,fontsize=17)
fig.savefig(out/'analytic-annulus.png',dpi=160);plt.close(fig)

# All mesh lines are loaded from the actual accepted solver inputs.
specs=[('原困难 Cut-cell',root/'build/laminar-fixtures/annulus-l7.solver.cm2d'),('GEOS 独立切割参考',study/'geos96.solver.cm2d'),('自研混合研究网格',study/'preserved-hybrid.solver.cm2d'),('cfMesh 完整参考（物面有偏移）',study/'cfmesh-native.solver.cm2d')]
fig,axes=plt.subplots(2,4,figsize=(18.5,9),layout='constrained')
for col,(name,path) in enumerate(specs):
 m=nm.read_cm2d(path);polys=[[m.vertices[v] for v in c.vertices] for c in m.cells]
 colors=['#a8d9be' if col==2 and max(math.hypot(*v) for v in p)<.522 else '#dceaf4' for p in polys]
 boundaries=[(m.vertices[e.v0],m.vertices[e.v1]) for e in m.edges if e.neighbour<0]
 for row in (0,1):
  ax=axes[row,col];ax.set_facecolor('#f2f2ef');ax.add_collection(PolyCollection(polys,facecolors=colors,edgecolors='#526d7d',linewidths=.3 if row==0 else .65));ax.add_collection(LineCollection(boundaries,colors='#263a45',linewidths=1.))
  ax.set_aspect('equal');ax.set_xlabel('x / m');ax.set_ylabel('y / m')
  if row==0:ax.set_xlim(-1.03,1.03);ax.set_ylim(-1.03,1.03);ax.set_title(f'{name}\n{len(m.cells):,} 格',fontproperties=font,fontsize=14)
  else:ax.set_xlim(-.55,-.39);ax.set_ylim(.15,.31);ax.set_title('同一内壁位置的放大',fontproperties=font,fontsize=12)
fig.suptitle('前三列保持同一物面；cfMesh 末列有几何近似；绿色为保留的近壁层',fontproperties=font,fontsize=16)
fig.savefig(out/'matched-meshes.png',dpi=160);plt.close(fig)

def read(path):
 with path.open() as f:return [{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
profiles=[('原困难网格',read(root/'outputs/cutcell-cross/original-coupled-strict.cells.csv'),read(study/'original-foam-coupled-v2/comparison.cells.csv')),
 ('混合研究网格',read(study/'preserved-hybrid-coupled.cells.csv'),read(study/'hybrid-external-coupled/comparison.cells.csv'))]
fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for col,(name,native,external) in enumerate(profiles):
 for rows,label,color in [(native,'自研 coupled','#2574ad'),(external,'成熟 coupled','#cc7b2a')]:
  gauge=model.errors(rows)['removedPressureGauge'];rad=np.array([math.hypot(x['x'],x['y']) for x in rows]);vel=np.array([(-x['y']*x['u']+x['x']*x['v'])/a for x,a in zip(rows,rad)]);pr=np.array([x['p']-gauge for x in rows])
  axes[0,col].scatter(rad,vel,s=2,alpha=.13,color=color,label=label,rasterized=True)
  axes[1,col].scatter(rad,pr,s=2,alpha=.13,color=color,label=label,rasterized=True)
 axes[0,col].plot(r,v,color='#232d32',lw=1.7,label='真圆解析解');axes[1,col].plot(r,p,color='#232d32',lw=1.7,label='真圆解析解')
 axes[0,col].set_title(name,fontproperties=font,fontsize=15)
 for row in (0,1):axes[row,col].set_xlim(.495,1.005);axes[row,col].grid(alpha=.16);axes[row,col].set_xlabel('半径 r / m')
 axes[0,col].set_ylabel('切向速度 / (m/s)');axes[1,col].set_ylabel('压力 p/ρ（去常数）/ (m²/s²)')
 axes[1,col].set_yscale('symlog',linthresh=.06)
 axes[1,col].set_ylim(-3,2)
 axes[1,col].set_yticks([-1,-.1,0,.1,1],labels=['−1','−0.1','0','0.1','1'])
 axes[0,col].legend(markerscale=4,fontsize=10)
fig.suptitle('实际收敛场：压力用对称对数纵轴，同时显示主体分布与局部尖峰',fontproperties=font,fontsize=15)
fig.savefig(out/'converged-fields-versus-circle.png',dpi=170);plt.close(fig)
print('figures',out)
