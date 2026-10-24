import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

mu=3.5e-3
def k(R,L): return np.pi*R**4/(8*mu*L)
mm=1e-3
segs=[("g1 lobar",1.20,15.0),("g2 segmental",0.70,12.0),("g3 subseg",0.40,10.0)]
res=[(n,1.0/k(R*mm,L*mm)) for n,R,L in segs]
Rseg=sum(r for _,r in res)
R_bed=0.40*Rseg; R_up=0.05*Rseg; R_tot=R_up+Rseg+R_bed
phis=[r/R_tot for _,r in res]

fig,(ax1,ax2)=plt.subplots(1,2,figsize=(10,3.9))

# Panel A: flow ratio vs delta for each break location
delta=np.linspace(-0.95,2.0,400)
colors=['#1b7837','#7f7f7f','#c51b7d']
for (n,r),phi,c in zip(res,phis,colors):
    fr=R_tot/(R_tot - r + r/(1+delta))
    ax1.plot(delta*100,fr,color=c,lw=2,label=f"{n}  ($\\phi$={phi:.3f})")
ax1.axhline(1.0,color='k',lw=0.6,ls=':')
ax1.axhspan(0.9,1.1,color='k',alpha=0.06)
ax1.axvline(0,color='k',lw=0.6,ls=':')
ax1.set_xlabel(r"bridge conductance error  $\delta=(k_{br}-k_{ab})/k_{ab}$  (%)")
ax1.set_ylabel(r"$q_{tum}^{broken}/q_{tum}^{intact}$")
ax1.set_title("(a)  Delivered-flow sensitivity to bridge error")
ax1.legend(fontsize=8,loc='upper right')
ax1.set_ylim(0,2.0); ax1.set_xlim(-95,200)
ax1.text(-90,0.2,"bridge with\nbulk tissue\n($\\delta\\to-1$)",fontsize=7,color='#555')

# Panel B: phi by generation (bar) with tolerance annotation
names=[n for n,_ in res]
x=np.arange(len(names))
bars=ax2.bar(x,phis,color=colors,width=0.6)
ax2.set_xticks(x); ax2.set_xticklabels([n.split()[0] for n in names])
ax2.set_ylabel(r"resistance fraction  $\phi$")
ax2.set_title("(b)  Per-break tolerance regime")
ax2.set_ylim(0,0.72)
for xi,phi in zip(x,phis):
    dtol=0.10/phi*100
    ax2.text(xi,phi+0.02,f"$\\phi$={phi:.3f}\n$\\pm${dtol:.0f}% for\n10% err",
             ha='center',va='bottom',fontsize=7.5)
ax2.axhline(0.1,color='k',ls='--',lw=0.7)
ax2.text(2.35,0.115,"wide  |  narrow",fontsize=7,color='#555',ha='right')
ax2.text(-0.4,0.66,"proximal / large R",fontsize=7,color='#1b7837')
ax2.text(-0.4,0.61,"distal / small R",fontsize=7,color='#c51b7d')

plt.tight_layout()
plt.savefig("phi_figure.pdf",bbox_inches='tight')
print("figure written")
