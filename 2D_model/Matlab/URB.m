function [f,x,k,var] = URB(dx,var)

% Set up of variables
a1=var.hyp.a1;
a2=var.hyp.a2;
a3=var.hyp.a3;
gammap=var.hyp.gp;
gamman=var.hyp.gn;
gammacump=var.hyp.gcp;
gammacumn=var.hyp.gcn;
gammlimp=var.hyp.glp;
gammlimn=var.hyp.gln;

fs1=var.hyp.fs1;
ps1=var.hyp.ps1;
fs2=var.hyp.fs2;
ps2=var.hyp.ps2;
fs3=var.hyp.fs3;
ps3=var.hyp.ps3;
fm =var.hyp.fm;
pm =var.hyp.pm;
h  =var.hyp.h;

xa=var.est.x;       % Displacement at the previous step
fa=var.est.f;       % Force at the previous step
x=xa+dx;            % Displacement at the new step


%% Hyperelastic component

% Evolution of the degradation factors

if x >= 0
    gammacum=gammacump;
    if dx >= 0
        gammlimp=max(gammlimp,x/h);
        gamma=min(gammap+dx/(2*h),gammlimp);
    else
        gammlimp=max(gammlimp,xa/h);
        gamma=min(gammap-4*dx/h,gammlimp);
        gammacum=gammacum-dx/h;
    end
    var.hyp.gp = gamma;
    var.hyp.gcp= gammacum;
    var.hyp.glp=gammlimp;
else
    gammacum=gammacumn;
    if dx < 0
        gammlimn=max(gammlimn,abs(x/h));
        gamma=min(gammlimn,gamman-dx/(2*h));
    else
        gammlimn=max(gammlimn,xa/h);
        gamma=min(gammlimn,gamman+4*dx/h);
        gammacum=gammacum+dx/h;
    end
    var.hyp.gn = gamma;
    var.hyp.gcn= gammacum;
    var.hyp.gln=gammlimn;
end

%Degradation factors
ks1= exp(-fs1*gamma^ps1);
ks2= ks1*exp(-fs2*gamma^ps2);
ks3= ks1*exp(-fs3*gamma^ps3);
km = exp(-fm*gammacum^pm);

fhyp=ks1*km*a1*(x/h)-ks2*a2*(x/h)^3+ks3*a3*(x/h)^5; 

%% Hysteretic component

% Set up of variables
uy=var.hys.uy;
fy=var.hys.fy;
k=var.hys.k;
B=var.hys.B;
g=var.hys.gamma;
n=var.hys.n;
fmax=var.hys.phimax;
pphi=var.hys.pphi;
ui=var.hys.ui;
uip=var.hys.uip;
zo=var.hys.z;
phi=var.hys.phi;
upa=var.hys.up;
up=dx/h;
u=x/h;

z=zo;
zm=zo;
flag=0;
count=1;

while flag==0
    zp=(k/fy)*up*(1-(abs(zm)^n)*(B*sign(up*zm)+g-phi*sign(up)*(sign(zm)+sign(up))));
    zf=zo+zp;
    if abs(zf-z)<0.0000001
        flag=1;
    else
        z=zf;
    end
    zm=zo+0.5*zp;
    if count > 100   % To avoid nonconvergence
        zm=zo;      % use the previous value
%         disp(['No conv in step ' num2str(ppo)])
    end
    count=count+1;
end

R=fy*z;         % Hysteretic force at the new step

%% Evolution of the variables for temporary hardening
 
var1=abs(sign(up)-sign(upa));
var2=abs(sign(z)-sign(zo));

% Considering the previous change of direction
if var1==2
    uip=abs(u-uy-up);
end

if var2==2
    phi=fmax*(1-exp(-pphi*abs(uip/uy)));
end

%% Update of variables

var.hys.z=z;
var.hys.up=up;
var.hys.uy=uy;
var.hys.phi=phi;
var.hys.uip=uip;
var.hys.ui=ui;

%% Resultinf forces
fr=R+fhyp;
f=[fhyp;R;fr];
k=(fr-fa)/(x-xa);
var.est.x=x;
var.est.f=fr;

end

