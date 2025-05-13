function [ F, v, var ] = HDRB_3D( inc, var )
    
    du = [inc(1);inc(2)];
    dP = inc(3);
    
    [ F, var ] = TSM( dP, du, var);
    
    kv_i = var.cav.kvi;

    v = var.est.v;
    u = norm(var.est.u);
    
    a = var.Ri/var.Ro;
    R = sqrt(1+a^2)*var.Ro/4;
    kv=kv_i/(1+(var.fdv/pi^2)*((u/R)^2));


P = var.est.P;


if P <= 0
    v = P/kv;
else
    uc = var.cav.uc;
    um=max(v,var.cav.um);
    phi_m = var.cav.phim;
    a = var.cav.a;
    f_c = var.cav.fc;
    k= var.cav.k;
    tr=var.Tr;

    phi = phi_m*(1-exp(-a*((um-uc)/uc)));
    f_cn =f_c*(1-phi);
    u_cn=f_cn/kv;

    f_pc=var.cav.fpc;

    kd=(f_pc-f_cn)/(um-u_cn);
    
    if P <= f_cn
        v=P/kv;
        
    else if P <= f_pc
            
            v = (P-f_cn)/kd +u_cn;
        else
            v = uc-log(1-k*tr*(P/f_c-1))/k;
        end
    end
    
    var.cav.um=um;
    var.cav.fpc=max(f_pc,P);
    
    var.est.v=v;
    

end

