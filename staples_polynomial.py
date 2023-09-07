'''
this script calculate staples in polynomial. 
it use for generate coefficients for equation intersect of line and surface
'''

from  collections import defaultdict

class Poli:
    def __init__(self) -> None:
        self.p = defaultdict(list)

    def normalize(self):
        for key, item in self.p.items():
            nm = defaultdict(int)
            for a in item:
                k = a[0]
                s = a[1]
                s = ''.join(sorted(s))
                nm[s] += k
            item = []    
            for kk, ii in nm.items():    
                item.append((ii, kk))
            self.p[key] = item
    
    def format(self):            
        res = ""
        for key, item in self.p.items():
            s = f"float a{key} = "
            start = 0
            for a in item:
                if a[0] == 0:
                    continue
                c = f"{a[0]}."
                k = a[1]
                for e in k:
                    c = c + f"*{e}"
                if  a[0] > 0 and start > 0:
                    c = " + " + c
                s = s + c
                start +=1
            res = res + s + ";\n"    
        return res    
       
            



def addpo(po1, po2):
    res = Poli()
    for key, item in po1.p.items():
        res.p[key].extend(item)
    
    for key, item in po2.p.items():
        res.p[key].extend(item)    

    res.normalize()
    return res


def multipo(po1, po2):
    res = Poli()
    for key1, item1 in po1.p.items():
        for key2, item2 in po2.p.items():
            key = key1 + key2
            item = []
            for a1 in item1:
                for a2 in item2:
                    k = a1[0]*a2[0]
                    s = a1[1] + a2[1]
                    item.append((k, s))
            res.p[key].extend(item)
    res.normalize()        
    return res               

#clebsch surface
# k(x^3+y^3+z^3+t^3) - (x + y + z + t)^3
def clebsch_surface():
    x = Poli()
    x.p[0] = [(1, "a")]
    x.p[1] = [(1, "b")]
    y = Poli()
    y.p[0] = [(1, "c")]
    y.p[1] = [(1, "d")]
    z = Poli()
    z.p[0] = [(1, "e")]
    z.p[1] = [(1, "f")]
    t = Poli()
    t.p[0] = [(1, "t")]
    nega = Poli()
    nega.p[0] = [(-1, '')]
    k = Poli()
    k.p[0] = [(1, 'k')]

    r1 = addpo(addpo(addpo(x, y), z), t)  #(x+y+z+t)
    r1 = multipo(multipo(r1, r1), multipo(r1, nega))  #-(x+y+z+t)^3


    x = multipo(multipo(x, x), x)
    y = multipo(multipo(y, y), y)
    z = multipo(multipo(z, z), z)
    t = multipo(multipo(t, t), t)
    r2 = addpo(addpo(addpo(x, y), z), t)  #x^3+y^3+z^3
    r2 = multipo(r2, k) #k(x^3+y^3+z^3+t^3)
    
    r = addpo(r2, r1) #k(x^3+y^3+z^3+t^3) - (x + y + z + t)^3
    print(r.format()) 


#z(z^2 - t^2) - x(x^2 - 3y^2) = 0
#normal = -3x^2 + 3y^2, 6yx, 3z^2 - t^2
def surfase1():
    x = Poli()
    x.p[0] = [(1, "a")]
    x.p[1] = [(1, "b")]
    y = Poli()
    y.p[0] = [(1, "c")]
    y.p[1] = [(1, "d")]
    z = Poli()
    z.p[0] = [(1, "e")]
    z.p[1] = [(1, "f")]
    t = Poli()
    t.p[0] = [(1, "t")]
    nega = Poli()
    nega.p[0] = [(-1, '')]
    nega3 = Poli()
    nega3.p[0] = [(-3, '')]

    r1 = multipo(z, addpo(multipo(z, z), multipo(multipo(t, t), nega))) # z(z^2 - t^2)
    r2 = multipo(multipo(x, addpo(multipo(x, x), multipo(multipo(y, y), nega3))), nega) #- x(x^2 - 3y^2)
    r = addpo(r1, r2) #z(z^2 - t^2) - x(x^2 - 3y^2)
    print(r.format())


#(x+y+z-a)(xy+yz+zx) - kxyz = 0
#normal = (xy+yz+zx) + (y+z)(x+y+z-a) - kyz, (xy+yz+zx) + (x+z)(x+y+z-a) - kxz,  (xy+yz+zx) + (y+x)(x+y+z-a) - kxy
def surface2():
    x = Poli()
    x.p[0] = [(1, "a")]
    x.p[1] = [(1, "b")]
    y = Poli()
    y.p[0] = [(1, "c")]
    y.p[1] = [(1, "d")]
    z = Poli()
    z.p[0] = [(1, "e")]
    z.p[1] = [(1, "f")]
    t = Poli()
    t.p[0] = [(1, "t")]
    nega = Poli()
    nega.p[0] = [(-1, '')]
    k = Poli()
    k.p[0] = [(1, 'k')]

    r1 = addpo(addpo(addpo(x, y), z), multipo(t, nega)) #(x+y+z-a)
    r2 = addpo(addpo(multipo(x, y), multipo(y, z)), multipo(z,x)) #(xy+yz+zx)
    r3 = multipo(multipo(multipo(x, y), multipo(k, z)), nega) #kxyz
    r = addpo(multipo(r1, r2), r3) #(x+y+z-a)(xy+yz+zx) - kxyz = 0    
    print(r.format())

#clebsch_surface()
#surfase1()
surface2()






