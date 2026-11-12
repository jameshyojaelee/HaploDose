"""Independent centered-product correlation and a float32 forward-error enclosure.

Sequential-sum gamma bounds also bound pairwise summation. The enclosure allows
mean, subtraction, product, sum, sqrt and division rounding separately. This is
an arithmetic check, never an uncertainty interval for the scientific estimand.
"""
import math
import numpy as np


def verify():
    for y,expected in [([-1,0,1],1.),([1,0,-1],-1.),([-1,1,0],.5)]:
        actual=reference([-1,0,1],y)
        if abs(actual['value']-expected)>4*np.finfo(float).eps:raise AssertionError('Analytic Pearson fixture')
        if not actual['float32_lower']<=expected<=actual['float32_upper']:raise AssertionError('Arithmetic enclosure')


def reference(x,y):
    x,y=np.asarray(x,np.longdouble),np.asarray(y,np.longdouble)
    if x.ndim!=1 or x.shape!=y.shape or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Finite aligned vectors required')
    n=len(x);u=np.longdouble(np.finfo(np.float32).eps)/2
    def gamma(k):return k*u/(1-k*u)
    if (n+8)*u>=1:raise ValueError('Arithmetic bound unavailable at this length')
    mx,my=x.mean(),y.mean();dx,dy=x-mx,y-my
    sx,sy=(dx*dx).sum(),(dy*dy).sum();dot=(dx*dy).sum()
    if sx<=0 or sy<=0:raise ValueError('Constant vector has undefined Pearson correlation')
    value=dot/np.sqrt(sx*sy)
    ex=gamma(n+1)*np.abs(x).sum()/n+u*(np.abs(x)+abs(mx))
    ey=gamma(n+1)*np.abs(y).sum()/n+u*(np.abs(y)+abs(my))
    dxerr=np.sum(2*np.abs(dx)*ex+ex*ex)+gamma(n+1)*np.sum((np.abs(dx)+ex)**2)
    dyerr=np.sum(2*np.abs(dy)*ey+ey*ey)+gamma(n+1)*np.sum((np.abs(dy)+ey)**2)
    derr=np.sum(np.abs(dx)*ey+np.abs(dy)*ex+ex*ey)+gamma(n+1)*np.sum((np.abs(dx)+ex)*(np.abs(dy)+ey))
    # Extra operations cover norm rescaling used by SciPy (divide, sqrt,
    # multiplication, and final normalized sum) without assuming a BLAS order.
    if sx<=dxerr or sy<=dyerr:raise ArithmeticError('Float32 correlation numerically unidentified')
    lower=np.sqrt((sx-dxerr)*(sy-dyerr))*(1-gamma(8))
    upper=np.sqrt((sx+dxerr)*(sy+dyerr))*(1+gamma(8))
    candidates=[v/d for v in (dot-derr,dot+derr) for d in (lower,upper)]
    lo,hi=min(candidates),max(candidates)
    # A second independent binary64 sum implementation checks refinement.
    a,b=np.asarray(x,float),np.asarray(y,float)
    ac=a-math.fsum(a)/n;bc=b-math.fsum(b)/n
    second=math.fsum(ac*bc)/math.sqrt(math.fsum(ac*ac)*math.fsum(bc*bc))
    difference=abs(float(value)-second)
    budget=8*n*np.finfo(float).eps*(1+float(np.abs(dx*dy).sum()/np.sqrt(sx*sy)))
    if difference>budget:raise ArithmeticError('Independent double/extended correlation disagreement')
    return dict(value=float(value),float32_lower=float(lo),float32_upper=float(hi),
                double_extended_difference=difference,double_roundoff_budget=budget,n=n)
