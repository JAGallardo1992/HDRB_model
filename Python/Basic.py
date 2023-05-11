

def sign(val):
    if abs(val) ==0:
        res =0
    else:
        res = val/abs(val)
    return res


def rest_vec(a,b):
    na=len(a)
    nb=len(b)
    if na == nb:
        c=list()
        for i in range(na):
            c.append(a[i]-b[i])
    else:
        print('Las listas deben ser de la misma longitud')
        return
    return c

def sum_vec(d,e):
    nd=len(d)
    ne=len(e)
    if nd == ne:
        ff=list()
        for i in range(nd):
            ff.append(d[i]+e[i])
    else:
        print('Las listas deben ser de la misma longitud')
        return
    return ff

def pot_vec(q,n):
    w=list()
    for i in q:
        w.append(i**n)
    return w


def heaviside(value):
    if value > 0:
        result = 1
    else:
        result = 0
    return result

