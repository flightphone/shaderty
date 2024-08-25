#print(3.14159265359*2)
def cpoli(s):
    res = Poli()
    res.p[0] = [(1, s)]
    return res



def step2():
    a = cpoli("a")
    b = cpoli("b")
    c = cpoli("c")
    d = cpoli("d")
    e = cpoli("e")
    f = cpoli("f")
    #g = cpoli("g")
    #h = cpoli("h")
    #i = cpoli("i")
    #k = cpoli("k")
    #l = cpoli("l")
    #u = cpoli("u")
    v = cpoli("v")
    A = cpoli("A")  #g/a
    B = cpoli("B")  #i/d

    

    F = cpoli("F")
    C = cpoli("C")
    D = cpoli("D")
    G = cpoli("G")
    H = cpoli("H")

    K = cpoli("K") #K = G*G*a/d

    x = cpoli("x")
    y = cpoli("y")
    z = cpoli("z")
    up = C*x + D*y + F*z + G*v + H
    xp = a*up*up + b*up + c
    yp = d*v*v + e*v + f
    zp = xp - K*yp
    print(zp.format())

def step1():
    a = cpoli("a")
    b = cpoli("b")
    c = cpoli("c")
    d = cpoli("d")
    e = cpoli("e")
    f = cpoli("f")
    g = cpoli("g")
    h = cpoli("h")
    i = cpoli("i")
    k = cpoli("k")
    l = cpoli("l")
    u = cpoli("u")
    v = cpoli("v")
    A = cpoli("A")  #g/a
    B = cpoli("B")  #i/d

    x = a*u*u + b*u + c
    y = d*v*v + e*v + f
    z = g*u*u + h*u + i*v*v + k*v + l
    r1 = z - x*A - y*B
    print(r1.format())


def bezier(x, y, z):
    a = cpoli("a")
    b = cpoli("b")
    c = cpoli("c")
    d = cpoli("d")
    e = cpoli("e")
    f = cpoli("f")
    #g = cpoli("g")
    #h = cpoli("h")
    #i = cpoli("i")
    #k = cpoli("k")
    #l = cpoli("l")
    #u = cpoli("u")
    v = cpoli("v")
    A = cpoli("A")  #g/a
    B = cpoli("B")  #i/d

    

    F = cpoli("F")
    C = cpoli("C")
    D = cpoli("D")
    G = cpoli("G")
    H = cpoli("H")

    K = cpoli("K") #K = G*G*a/d
    t2 = Poli()
    t2.p[0] = [(2, '')]
    f1 = t2*C*G*a*x + t2*D*G*a*y + t2*F*G*a*z + t2*G*H*a + G*b-K*e
    f2 = C*C*a*x*x + t2*C*D*a*x*y + t2*C*F*a*x*z + t2*C*H*a*x + D*D*a*y*y + t2*D*F*a*y*z  + t2*D*H*a*y + F*F*a*z*z  + t2*F*H*a*z  + H*H*a + C*b*x + D*b*y + F*b*z  + H*b -K*f
    f3 = x - K*y - f2
    res =  d * f3 * f3 +  e*f3*f1 + f*f1*f1 - y*f1*f1
    return res


def bezier_surf():
    x, y, z = param2()
    res = bezier(x,y,z)
    print(res.format2())

    x, y, z = axis()
    res = bezier(x,y,z)
    print(res.normal())

#bezier_surf()
