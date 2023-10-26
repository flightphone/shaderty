'''
this script calculate staples in polynomial. 
it use for generate coefficients for equation intersect of line and surface
'''

from  collections import defaultdict
from collections import Counter
import sys
from poly import *

#clebsch surface
# k(x^3+y^3+z^3+t^3) - (x + y + z + t)^3
def clebsch_surface(x, y, z, t) -> Poli:
    k = Poli()
    k.p[0] = [(1, 'k')]
    u = x + y + z + t
    return k*(x*x*x + y*y*y + z*z*z + t*t*t) - u*u*u
    #k(x^3+y^3+z^3+t^3) - (x + y + z + t)^3


#z(z^2 - t^2) - x(x^2 - 3y^2) = 0
#normal = -3x^2 + 3y^2, 6yx, 3z^2 - t^2
def example_surfase(x, y, z, t) -> Poli:
    k = Poli()
    k.p[0] = [(3, '')]
    return z*(z*z - t*t) - x*(x*x - k*y*y)
    #z(z^2 - t^2) - x(x^2 - 3y^2)

#CAYLEY SURAFCE https://mathcurve.com/surfaces.gb/cayley/cayley.shtml
#(x+y+z-a)(xy+yz+zx) - kxyz = 0
#normal = (xy+yz+zx) + (y+z)(x+y+z-a) - kyz, (xy+yz+zx) + (x+z)(x+y+z-a) - kxz,  (xy+yz+zx) + (y+x)(x+y+z-a) - kxy
def gayley(x, y, z, t) -> Poli:
    k = Poli()
    k.p[0] = [(1, 'k')]
    return (x+y+z-t)*(x*y + y*z + z*x) - k*x*y*z
    #(x+y+z-a)(xy+yz+zx) - kxyz = 0
    

#MÖBIUS SURFACE
#https://mathcurve.com/surfaces.gb/mobiussurface/mobiussurface.shtml
def mobius(x, y, z, t) -> Poli:
    k2 = Poli()
    k2.p[0] = [(2, '')]
    yz = y - z 
    xt = x + t 
    r = y*(yz*yz - xt*xt) + yz * xt* x * k2
    return r
    #y((y-z)^2 - (x+a)^2) + 2x(y-z)(x+a) = 0


#MONKEY SADDLE
#a^2 z = x (x^2 - 3y^2)  //https://mathcurve.com/surfaces.gb/selle/selle.shtml
def monkey(x, y, z, t) -> Poli:
    k = Poli()
    k.p[0] = [(3, '')]
    return x*(x*x - k*y*y) - t*t*z
    #x (x^2 - 3y^2) - a^2 z
    
#ROMAN SURFACE https://mathcurve.com/surfaces.gb/romaine/romaine.shtml
def roman(x, y, z, t) ->Poli:
    k = Poli()
    k.p[0] = [(2, '')]
    return x*x*y*y + y*y*z*z + x*x*z*z - k*t*x*y*z
    #x^2 y^2 + y^2 z^2 + x^2 z^2 - 2axyz


#GOURSAT SURFACE https://mathcurve.com/surfaces.gb/goursat/goursat.shtml
def goursat(x, y, z, t) ->Poli:
    k = Poli()
    k.p[0] = [(2, '')]
    return k*(x*x*y*y + y*y*z*z + x*x*z*z) - t*t*(x*x + y*y + z*z) - t*t*t*t
    #2(x^2 y^2 + y^2 z^2 + x^2 z^2) - t^2(x^2 + y^2 + z^2) - t^4

#DUPIN CYCLIDE https://mathcurve.com/surfaces.gb/cycliddedupin/cyclidededupin.shtml
def cyclide(x, y, z, t) ->Poli:
    k = Poli()
    k.p[0] = [(4, '')]

    u = Poli()
    u.p[0] = [(1, 'u')]

    w = Poli()
    w.p[0] = [(1, 'w')]

    v = Poli()
    v.p[0] = [(1, 'v')]

    return (x*x + y*y + z*z + u*u - w*w)*(x*x + y*y + z*z + u*u - w*w) - k*((t*x - v*w)*(t*x - v*w) + u*u*y*y)
    #(x^2 + y^2 + z^2 + u^2 - w^2)^2 - 4((tx - v*w)^2 + u^2*y^2)
def test(x, y, z, t):
    return (x*x - x)*(y*y - y) - z
    
#clebsch_surface()
#example_surfase()
#gayley()
#mobius()
#monkey()
#roman

sys.stdout = open("code.txt", "w")
polynom = test
x, y, z, t = param()
r = polynom(x, y, z, t)
print(r.format())

x, y, z = axis()
r = polynom(x, y, z, t)
print(r.normal())






