import math, time
E=1.602176634e-19; KB=1.380649e-23; NA=6.02214076e23; EPS0=8.8541878188e-12
T=298.15; EPSR=78.30; LAM=100e-9; GB=14.0
def logW(c, psi_mV, a_nm, A_zJ, retard=True, hydro=True, N=1<<14, umin=1e-10, umax=1e6):
    kT=KB*T; a=a_nm*1e-9; A=A_zJ*1e-21; psi=psi_mV*1e-3
    kap=math.sqrt(2*NA*1000*c*E*E/(EPSR*EPS0*kT))
    s0,s1=math.log(umin),math.log(umax); ds=(s1-s0)/N
    ln_num=[];ln_den=[]
    for i in range(N+1):
        s=s0+i*ds; u=math.exp(s); h=u*a
        va=-A*a/(12*h)/kT
        if retard: va/= (1+GB*h/LAM)
        vr=2*math.pi*EPSR*EPS0*a*psi*psi*math.log1p(math.exp(-kap*h))/kT
        b=(6*u*u+13*u+2)/(6*u*u+4*u) if hydro else 1.0
        w=(1 if i in (0,N) else (4 if i%2 else 2))
        base=math.log(w*b*u)-2*math.log(2+u)   # du = u ds
        ln_num.append(base+va+vr); ln_den.append(base+va)
    def lse(v):
        m=max(v); return m+math.log(sum(math.exp(x-m) for x in v))
    tail=math.log(1/(2+umax)) - math.log(ds/3)
    num=lse(ln_num+[tail]); den=lse(ln_den+[tail])
    return (num-den)/math.log(10)
if __name__=="__main__":
    for args in [(0.001,30,200,10),(0.01,30,200,10),(0.05,30,200,10),(0.1,30,200,10),(0.3,30,200,10),(0.01,15,100,10),(0.001,40,500,5),(0.05,20,50,20)]:
        t=time.time(); r=[logW(*args,N=n) for n in (1<<13,1<<14,1<<15)]
        print(args,[round(x,6) for x in r],'noRet',round(logW(*args,retard=False),4),'noHyd',round(logW(*args,hydro=False),4),'none',round(logW(*args,retard=False,hydro=False),4),f'{time.time()-t:.2f}s')
