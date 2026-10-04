"""Analytic flow solutions used to draw reference curves."""
import math

def sod(x,time,gamma=1.4,split=.5):
    """Exact self-similar Sod Riemann state, independently bracketed star pressure."""
    left=(1.,0.,1.);right=(.125,0.,.1)
    def wave(p,state):
        rho,u,pk=state;a=math.sqrt(gamma*pk/rho)
        if p>pk:return (p-pk)*math.sqrt((2/((gamma+1)*rho))/(p+(gamma-1)/(gamma+1)*pk))
        return 2*a/(gamma-1)*((p/pk)**((gamma-1)/(2*gamma))-1)
    low,high=.1,1.
    for _ in range(100):
        middle=(low+high)/2
        if wave(middle,left)+wave(middle,right)>0:high=middle
        else:low=middle
    ps=(low+high)/2;us=(wave(ps,right)-wave(ps,left))/2
    if time==0:return (1.,0.,1.) if x<split else (.125,0.,.1)
    xi=(x-split)/time;al=math.sqrt(gamma);ar=math.sqrt(gamma*.1/.125);astar=al*ps**((gamma-1)/(2*gamma))
    if xi<=-al:return left
    if xi<=us-astar:
        u=2/(gamma+1)*(al+xi);a=2/(gamma+1)*(al-(gamma-1)*xi/2)
        return (a/al)**(2/(gamma-1)),u,(a/al)**(2*gamma/(gamma-1))
    if xi<=us:return ps**(1/gamma),us,ps
    shock=ar*math.sqrt((gamma+1)/(2*gamma)*ps/.1+(gamma-1)/(2*gamma))
    if xi<shock:
        ratio=ps/.1;k=(gamma-1)/(gamma+1)
        return .125*(ratio+k)/(k*ratio+1),us,ps
    return right


def linear_euler_fourier_mode(time, conductivity, gamma=1.4, gas_r=1., rho=1., temperature=1., wavelength=1., amplitude=1e-5, viscosity=0.):
    """Continuum linearized Euler-Fourier mode, independent of spatial stencils.

    drho=A*cos(wx), u=B*sin(wx), dT=C*cos(wx). Initially isobaric
    to first order: A=-rho*amplitude, B=0, C=T*amplitude.
    Evaluate exp(time*M) by scaling/squaring of its convergent Taylor series.
    The nonlinear PDE differs by O(amplitude**2); this is a small-signal test.
    """
    w=2*math.pi/wavelength
    alpha=conductivity/(rho*gas_r/(gamma-1))
    m=[[0.,-rho*w,0.], [gas_r*temperature*w/rho,-4/3*viscosity/rho*w*w,gas_r*w],
       [0.,-(gamma-1)*temperature*w,-alpha*w*w]]
    norm=max(sum(abs(x*time) for x in row) for row in m)
    squarings=max(0,math.ceil(math.log2(norm/.5))) if norm else 0
    h=time/(2**squarings)
    a=[[x*h for x in row] for row in m]
    identity=lambda:[[float(i==j) for j in range(3)] for i in range(3)]
    def multiply(left,right):
        return [[math.fsum(left[i][k]*right[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    result=identity();term=identity()
    for n in range(1,64):
        term=[[x/n for x in row] for row in multiply(term,a)]
        result=[[x+y for x,y in zip(row,delta)] for row,delta in zip(result,term)]
        if max(abs(x) for row in term for x in row)<math.ulp(1.)/8:break
    else:raise ValueError('matrix exponential did not converge')
    for _ in range(squarings):result=multiply(result,result)
    initial=(-rho*amplitude,0.,temperature*amplitude)
    return tuple(math.fsum(x*y for x,y in zip(row,initial)) for row in result)

