function [ F, var, iter, s] = TSM( dP, du, var)

u_a = var.est.u;
P_a = var.est.P;

P = P_a+dP;
u = u_a+du;

var.est.P = P;

tol = var.tol;
max_iter =var.max_iter;

h=var.hyp.h;
th=var.TSM.th;
s_a=var.TSM.s;
kt=var.TSM.kt;
F=var.TSM.F;
fs=var.est.f;

Xn=[th;s_a;fs;F];
error=1;
iter=0;
while error > tol
    s=u-h*th;
    inc=s-s_a;
    [fs,var_n] = BIHDRB3( inc,var );       
    F=fs-P*th;
    th=real((F*h+P*s)/(kt-P*h));
    Xt=[th;s;fs;F];
    e2=sum((Xt-Xn).^2);
    e1=heaviside(max_iter-iter);
    error=min(e1,e2);
    Xn=Xt;
    iter=iter+1;
end

var=var_n;
var.est.u=u;
var.TSM.th=th;
var.TSM.s=s;
var.TSM.F=F;
var.est.f=fs;
end

