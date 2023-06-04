#______________________________Gravipy for tensors and metric___________________________

# #________Basics gravipy____________

from gravipy.tensorial import * # import GraviPy package
from sympy import init_printing
import inspect
init_printing()

# # define some symbolic variables
t, r, theta, phi, M = symbols('t, r, \\theta, \phi, M')  #real symbols
# # create a coordinate four-vector object instantiating 
# # the Coordinates class
x = Coordinates('\chi', [t, r, theta, phi])

# # define a matrix of a metric tensor components _Schwarschild
Metric = diag(-(1-2*M/r), 1/(1-2*M/r), r**2, r**2*sin(theta)**2)  #metric ??
# # create a metric tensor object instantiating the MetricTensor class
g = MetricTensor('g', x, Metric)

#print(x(-All))
# #print(g(All,All))

# #_______Predefined tensor classes________
# print([cls.__name__ for cls in vars()['Tensor'].__subclasses__()])

Ga = Christoffel('Ga', g)
# Ga(1, 2, 1)

# print(Ga.components)
# help(Christoffel)

Ri = Ricci('Ri', g)
# Ri(All, All)
Rm = Riemann('Rm', g)

from IPython.display import display, Math
from sympy import latex
for i, j, k, l in list(variations(range(1, 5), 4, True)):
    if Rm(i, j, k, l) != 0 and k<l and i<j:
        display(Math('R_{'+str(i)+str(j)+str(k)+str(l)+'} = '+ latex(Rm(i, j, k, l))))

# #______Geodesics_________

tau = Symbol('\\tau')
w = Geodesic('w', g, tau)
w(All).transpose()

Parametrization.info()
Parametrization.deactivate(x)
Parametrization.info()
T = Tensor('T', 2, g)
T.partialD(1, 2, 1, 3)
print(T.partial_derivative_components)

C = Tensor('C', 0, g)

for k in C.covariant_derivative_components:
    display(Math(str(k) + ': '  + latex(C.covariant_derivative_components[k])))
    
not any([g.covariantD(i, j, k).simplify() for i, j, k in list(variations(range(1, 5), 3, True))])






