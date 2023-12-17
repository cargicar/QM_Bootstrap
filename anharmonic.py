import numpy as np

def xs(s,E,g,x2):
    ''' Equation 6 Hartnoll.
    Input: x2 = x**2
    Returns x**s
    Notice that this can be improved with dynamical programming'''
    if s%2==1 or s < 0:
        return 0
    elif s==0:
        return 1
    elif s==2:
        return x2
    return (4*(s-3)*E*xs(s-4,E,g,x2)+(s-3)*(s-4)*(s-5)*xs(s-6,E,g,x2)-4*(s-2)*xs(s-2,E,g,x2))/(4*g*(s-1))

xsarray=[1,0]
def xs_dyn(s,E,g,x2):
    ''' Equation 6 Hartnoll.
    Input: x2 = x**2
    Returns x**s
   Dynamic programming approach '''
    if s==2:
        xsarray.append(x2)
        return xsarray
    else:
        xsarray.append((4*(s-3)*E*xs(s-4,E,g,x2)+(s-3)*(s-4)*(s-5)*xs(s-6,E,g,x2)-4*(s-2)*xs(s-2,E,g,x2))/(4*g*(s-1)))
    return xsarray

def M(size,E,g,x2):
    m=np.zeros((size, size))
    for i in range(size):
        for j in range(size):
            m[i,j]=xs(i+j,E,g,x2)
    return m

def M_dyn(size,E,g,x2):
    '''under construction'''
    m=np.zeros((size, size))
    for i in range(size):
        for j in range(size):
            m[i,j]=xs_dyn(i+j,E,g,x2)[-1]
    return m

def pert_e(w, g):
    a0=  0.667986259
    a1=  0.143668783
    a2= -0.008627565
    return (g/2)**(1/3)*(a0+a1*(g/(2*w**3))**(-2/3)+a2*(g/(2*w**3))**(-4/3) )
# np.linalg.cholesky(A)

def e_ben(l):
    ans=[3/4, -21/8, 333/16, -30885/128, 916731/256, -65518401/1024]
    ls=[(l/2)**i*ans[i] for i in range(len(ans))]
    return 1/2+sum(ls)
