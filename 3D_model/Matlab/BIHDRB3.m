function [ F, var2 ] = BIHDRB3( inc,var2 )

uy = var2.hys.uy;
fy = var2.hys.fy;
k = var2.hys.k;
B = var2.hys.B;
gamma = var2.hys.gamma;
n = var2.hys.n;

pphi = var2.hys.pphi;
uip = var2.hys.uip;
zo = var2.hys.z;
phi = var2.hys.phi;
fmax = var2.hys.phimax;

up = inc/var2.hyp.h;
ua = var2.hys.u;
u = ua+inc;

z=zo;
zm=zo;
flag=0;
count=1;

while flag == 0
    if (zm(1)==0) && (zm(2)==0)
        zm(1)=1e-50;
        zm(2)=1e-50;
    end
    zp_x = (k/fy)*up(1)-(k/fy)*zm(1)*((B-phi)*abs(up(1)*zm(1))+(gamma-phi)*up(1)*zm(1)+(B-phi)*abs(up(2)*zm(2))+(gamma-phi)*up(2)*zm(2))*(zm(1)^2+zm(2)^2)^(n/2-1);
    zp_y = (k/fy)*up(2)-(k/fy)*zm(2)*((B-phi)*abs(up(1)*zm(1))+(gamma-phi)*up(1)*zm(1)+(B-phi)*abs(up(2)*zm(2))+(gamma-phi)*up(2)*zm(2))*(zm(1)^2+zm(2)^2)^(n/2-1);
   
    zp = [zp_x;zp_y];
    zf =zo+zp;
   
    if norm(zf-z)<0.0000001
        flag=1;
    else
        z=zf;
    end
    zm=zo+0.5*zp;
    if count > 100   % To avoid nonconvergence
        zm=zo;      % use the previous value
    end
    count=count+1;
end

F1 = fy*z;         % Actual hyteretic force

%% Evolution of the variables for isotropic hardening
N_u = norm(u);
N_ua = norm(ua);

if N_u < N_ua
    uip=max(N_ua-uy,uip);
    phi=fmax*(1-exp(-pphi*abs(k*uip/(var2.hyp.h*fy))));
end


a1=var2.hyp.a1;
a2=var2.hyp.a2;
a3=var2.hyp.a3;

fs1=var2.hyp.fs1;
ps1=var2.hyp.ps1;
fs2=var2.hyp.fs2;
ps2=var2.hyp.ps2;
fs3=var2.hyp.fs3;
ps3=var2.hyp.ps3;
fm =var2.hyp.fm;
pm =var2.hyp.pm;
h  =var2.hyp.h;

gu_n=norm(u)/h;
gu_a=norm(ua)/h;

um = norm(u);
if um==0
    u_n=[0;0];
else
    u_n=u/um;
end

dir=var2.hyp.dir;



ip_a=ua'*dir;
ip_n=u'*dir;
   
[~,ui_an]=find(ip_a==max(ip_a));
[~,ui_n]=find(ip_n==max(ip_n));

if length(ui_n)==1
    ui_a=ui_n;
elseif isempty(ui_n) ~=1
        ui_a=length(ui_n)+floor(length(ui_n)/2)-floor(length(ui_n)/2)*2;
else
        ui_a=length(ui_an)+floor(length(ui_an)/2)-floor(length(ui_an)/2)*2;
end


gammalim = var2.hyp.gammal(ui_a);
gamma = var2.hyp.gamma(ui_a);
gammac = var2.hyp.gammac(ui_a);


gammalimh = var2.hyp.gammalh(ui_a);


%% Evolution of the degradation factors
vec_gl = var2.hyp.gammal;
vec_g  = var2.hyp.gamma;
vec_gc = var2.hyp.gammac;

dg=gu_n-gu_a;

gammalim=max(gammalim,gu_n);
gamma_trial=gamma+(9*sign(dg)-7)*(dg/4);

gamma_a=gamma;

gamma=min(gammalim,gamma_trial);
gammac=gammac-(1-sign(dg))*(dg/2);

gammalimh = max(gammalimh,gu_n);

d_gamma=gamma-gamma_a;

vec_gl(ui_a) = gammalim;
vec_gc(ui_a) = gammac;

ui=dir(:,ui_a);
ang=acos(ui'*dir);
fcpl=var2.hyp.fc;
gl_trial  = (fcpl+(1-fcpl).*cos(ang)).*gammalimh;
gc_trial = (fcpl+(1-fcpl).*cos(ang)).*gammac;
dg_trial  = (fcpl+(1-fcpl).*cos(ang)).*d_gamma;


var2.hyp.gammal = max(gl_trial,vec_gl);
var2.hyp.gamma = min(var2.hyp.gammal,vec_g+dg_trial);
var2.hyp.gammac = max(gc_trial,vec_gc);



var2.hyp.gammalh(ui_a) = gammalimh;


gamma = var2.hyp.gamma(ui_a);
gammac = var2.hyp.gammac(ui_a);

ks1= real(exp(-fs1*gamma^ps1));
ks2= real(ks1*exp(-fs2*gamma^ps2));
ks3= real(ks1*exp(-fs3*gamma^ps3));

km = exp(-fm*gammac^pm);
%% Hyperelastic component

F2=(a1*km*ks1)*(um*u_n)-(a2*ks2)*(um^3*u_n)+(a3*ks3)*(um^5*u_n);
F=F1+F2;
%% Update of variables
var2.hys.u=u;
var2.hys.z=z;
var2.hys.up=inc;
var2.hys.uy=uy;
var2.hys.phi=phi;
var2.hys.uip=uip;

end

