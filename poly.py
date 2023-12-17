'''
this script calculate staples in polynomial. 
it use for generate coefficients for equation intersect of line and surface
'''

from  collections import defaultdict
from collections import Counter


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
    
    def df(self, v):
        res = "0."
        for item in self.p[0]:
            k = item[0]
            if k == 0:
                continue
            s = item[1]
            c = Counter(s)
            deg = c.get(v, 0)
            if deg == 0:
                continue
            k *= deg
            c[v] -= 1
            val = f"{k}."
            if k > 0:
                val = "+" + val
            for leter, count in c.items():
                if leter in {"x", "y", "z"}:
                    leter = "pos." + leter
                for i in range(count):
                    val = val + "*" + leter
            res = res + val
        return res            




    def normal(self):
        dx = self.df("x")
        dy = self.df("y")
        dz = self.df("z")
        res = f"nor = vec3({dx}, {dy}, {dz});"
        return res

    
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
       
    def __add__(self, po2):        
        return addpo(self, po2)
    
    def __sub__(self, po2):
        po1 = self
        n = Poli()
        n.p[0] = [(-1, '')]   
        return addpo(self, multipo(po2, n))
    
    def __neg__(self):
        n = Poli()
        n.p[0] = [(-1, '')]   
        return multipo(self, n)
    
    def __mul__(self, other):
        return multipo(self, other)




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

def axis():
    x = Poli()
    x.p[0] = [(1, "x")]
    y = Poli()
    y.p[0] = [(1, "y")]
    z = Poli()
    z.p[0] = [(1, "z")]
    return x, y, z

def param():
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
    return x, y, z, t


a = [255, 74, 45]
for i, e in enumerate(a):
    a[i] = e/255
print (a)    