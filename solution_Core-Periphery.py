from scipy import optimize
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def Y1(x):
    out = mu*f*x[0] + (1.-mu)/2.
    return out

def Y2(x):
    out = mu*(1.-f)*x[1] + (1.-mu)/2.
    return out

def G1(x):
    out = ( f*x[0]**(1-sigma) + (1-f)*(x[1]*T )**(1-sigma) )**(1/(1-sigma))
    return out

def G2(x):
    out =( f*( x[0]*T )**(1-sigma) + (1-f)*x[1]**(1-sigma) )**(1/(1-sigma))
    return out

def w1(x): # as defined in Fujita et al
    out =  ( Y1(x)*G1(x)**(sigma-1) + Y2(x)*G2(x)**(sigma-1)*(T**(1-sigma)) )**(1/sigma)
    return out

def w2(x): # as defined in Fujita et al
    out = ( Y1(x)*G1(x)**(sigma-1)*(T**(1-sigma)) + Y2(x)*G2(x)**(sigma-1) )**(1/sigma)
    return out

def fun(x):
    return ( ( x[0] -w1(x) )**2 +( x[1] -w2(x) )**2 )

def func(x):
    out = [x[0]-w1(x)]
    out.append( x[1]-w2(x) )
    return out

sigma=4.
T=1.3
tau=1./T
mu=0.4


x=[1,1]

min=x
minimum=0.1

npoints=500

fig = plt.figure()
fig.subplots_adjust(hspace=0.4, wspace=0.4)
for i in range(1, 2):
    T=T+0.*i
    ax = fig.add_subplot(1, 1, i)
        #    ax.text(0.5, 0.5, str((2, 3, i)),
        #    fontsize=18, ha='center')

    a=[]
    b=[]
    for ii in range(41):
        f=ii/40.
        minimum=0.1
        min=[1,1]
    
        x02 = fsolve(func,x)
        a.append(f)
        b.append(x02[0]*G1(x02)**(-mu)-x02[1]*G2(x02)**(-mu))


    plt.axhline(y=0.,linewidth=1, color='k')
    plt.axvline(x=0.5,linewidth=1, color='k')

    plt.ylabel(r'$w_1-w_2$',fontsize=20)
    plt.xlabel(r'$f=L_1/\mu$',fontsize=20)

    plt.plot(a,b,'b-',lw=5,alpha=0.5)
    plt.xlim(0, 1)

#plt.plot(xx1,yy1,'b-',lw=8,alpha=0.5)
plt.savefig('lowT_costs_sigma=4_mu=0.3.pdf', format='pdf',bbox_inches='tight')
#plt.tight_layout()
#plt.show()

