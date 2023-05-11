def URB(x,var):
    
    import math
    from Basic import sign
    
    # Set up of variables
    a1 = var.get('a1',0) 
    a2 = var.get('a2',0)
    a3 = var.get('a3',0)
    gp = var.get('gp',0.0)
    gn = var.get('gn',0.0)
    gcp = var.get('gcp',0)
    gcn = var.get('gcn',0)
    glp = var.get('glp',0)
    gln = var.get('gln',0)
    
    fs1 = var.get('fs1',0)
    ps1 = var.get('ps1',0)
    fs2 = var.get('fs2',0)
    ps2 = var.get('ps2',0)
    fs3 = var.get('fs3',0)
    ps3 = var.get('ps3',0)
    fm = var.get('fm',0)
    pm = var.get('pm',0)
    
    uy = var.get('uy',0) 
    fy = var.get('fy',0) 
    B = var.get('B',0) 
    n = var.get('n',0) 
    fmax = var.get('fmax',0) 
    pphi = var.get('pphi',0) 
    ui = var.get('ui',0) 
    uip = var.get('uip',0) 
    zo = var.get('z',0) 
    phi = var.get('phi',0) 
    upa = var.get('up',0) 

    try:
        h=var['h']
        if h <= 0:
            print('The height must be larger than 0')
            return 'Failed',0,0
    except:
        print('Please enter the height of the rubber')
        return 'Failed',0,0
    
    
    g = 1.0-B
    k = var.get('k',fy/uy)
    xa = var.get('x',0.0)  
    dx = x-xa

    # Evolution of the degradation factors

    if (x >= 0.0):
        gc = gcp
        if dx >= 0.0:
            glp = max(glp,x/h)
            gm = min(gp+dx/(2.0*h),glp)
        else:
            glp = max(glp,xa/h)
            gm = min(gp-4.0*dx/h,glp)
            gc =gc-dx/h;
        var['gp'] = gm
        var['gcp'] = gc
        var['glp'] = glp
    else:
        gc=gcn
        if (dx < 0):
            gln = max(gln,abs(x/h))
            gm = min(gln,gn-dx/(2.0*h))
        else:
            gln = max(gln,xa/h)
            gm = min(gln,gn + 4.0 *dx/h)
            gc=gc+dx/h
        var['gn'] = gm
        var['gcn'] = gc
        var['gln'] = gln
    

    # Degradation factors
    ks1 = math.exp(-fs1*gm**ps1)
    ks2 = ks1*math.exp(-fs2*gm**ps2)
    ks3 = ks1*math.exp(-fs3*gm**ps3)
    km = math.exp(-fm*gc**pm)

    # Hyperelastic component

    fhyp=ks1*km*a1*(x/h)-ks2*a2*(x/h)**3.0+ks3*a3*(x/h)**5.0

    # Hysteretic component

    up = dx/h
    u = x/h
    z = zo
    zm = zo
    tol=10**(-7)
    count = 0
    flag = 1

    while (flag ==1) :
        zp=(k/fy)*up*(1-(abs(zm)**n)*(B*sign(up*zm)+g-phi*sign(up)*(sign(zm)+sign(up))))
        zf=zo+zp
        if (abs(zf-z) < tol):
            flag = 0
        else:
            z=zf
        zm=zo+0.5*zp
        if (count > 9):
            zm=zo      # use the previous value
        count += 1


    R = fy*z         # Actual hyteretic force
    Ro = fy*zo       # Previous hysteretic force


    # Evolution of the variables for temporary hardening
    var1=abs(sign(up)-sign(upa))
    var2=abs(sign(z)-sign(zo))

    # Considering the maximum historical deformation
    # for this case, divide pphi/2
    # if (var1 == 2):
    #     uipn = abs(u-uy-up)
    #     uip = max(uipn,uip)
    #     phin = fmax*(1.0-math.exp(-pphi*abs(k*uip/fy)))
    #     phi = max(phi,phin)
    #     uy = u-R/k+fy/(k*(1.0-4.0*phi))*sign(up)
    #     ui = u-up
    # if (var2 == 2):
    #     uy = u+R/k-fy/(k*(1.0-4.0*phi))
    #     ui = u-R/(R-Ro)*up
    
    if (var1==2):
        uip=abs(u-uy-up)

    if (var2==2):
        phi=fmax*(1-math.exp(-pphi*abs(uip/uy)))

    frint = fhyp+R
    
    # Update of variables
    var['k'] = k
    var['z'] = z
    var['up'] = up
    var['uy'] = uy
    var['phi'] = phi
    var['uip'] = uip
    var['ui'] = ui
    var['x'] = x
    var['f'] = frint
    return frint , var