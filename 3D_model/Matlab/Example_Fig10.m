%% Valos a probar el modelo de Bouc Wen
clear 
clc
close all

%% Input

Info=load('Data/Fig10.txt'); % Change load pattern
var.fdv = 5.1;  % axial stiffness factor

%% Parameters

n_dir = 360;
dt=0.01;
A=323976.74240144744/1000;
ux=Info(:,1);
uy=Info(:,2);
P=Info(:,3);

var.hys.k=1.8278;
var.hys.fy=52.9144;
var.hys.uy=var.hys.fy/var.hys.k;
var.hys.B=0.90;
var.hys.gamma=1-var.hys.B;
var.hys.n=0.3026;
var.hys.phimax=0.1147;
var.hys.pphi=2*0.0101;
var.hys.gamma=1-var.hys.B;


var.hyp.a1= 0.757116079670000000;
var.hyp.a2= 0.011748698800000;
var.hyp.a3= 1.644804500000000e-04;
var.hyp.fs1= 0.262220930000000;
var.hyp.ps1= 0.510406980000000;
var.hyp.fs2= 10.049106000000000;
var.hyp.ps2= 0.295731500000000;
var.hyp.fs3= 14.733708000000000;
var.hyp.ps3= 0.104399000000;
var.hyp.fm = 0.153779286000000;
var.hyp.pm = 0.240096460000000;
var.hyp.fc = 0.5;
var.hyp.h  =204;
var.h  =204;
var.hys.k=var.hyp.h*var.hys.fy/var.hys.uy;

var.Ri = 100;
var.Ro = 650;
var.Tr = 6;  
var.cav.kv = 850; 
var.cav.kvi = 850;
var.TSM.kt=51817730.946643;

var.cav.k = 20; 
var.cav.phim = 0.75; 
var.cav.a = 1.0;

%%  Update variables
gammalim=1.5;
gamma=1.5;

index=1:n_dir;
index=[index index index];

angle=(0:360/n_dir:360*(n_dir-1)/n_dir)*pi/180;
var.hyp.dir=[cos(angle);sin(angle)];
var.hyp.gamma =zeros(1,n_dir);
var.hyp.gammal =zeros(1,n_dir);
var.hyp.gammac =zeros(1,n_dir);
dir=var.hyp.dir;

ui=dir(:,1);
n_cent = n_dir+1;
var.hyp.gammal(1)=gammalim;
var.hyp.gamma(1)=gamma;

for i=1:(n_dir-(2-rem(n_dir,2)))/2
    uip = dir(:,index(n_cent+i));
    angp = acos((ui'*uip)/(norm(ui)*norm(uip)));
    var.hyp.gammal(index(n_cent+i)) = (0.5+cos(angp)/2)*gammalim;
    var.hyp.gamma(index(n_cent+i))  = (0.5+cos(angp)/2)*gamma;
    
    uin = dir(:,index(n_cent-i));
    angn = acos((ui'*uin)/(norm(ui)*norm(uin)));
    var.hyp.gammal(index(n_cent-i)) = (0.5+cos(angn)/2)*gammalim;
    var.hyp.gamma(index(n_cent-i))  = (0.5+cos(angn)/2)*gamma;
end

ui=dir(:,n_dir/2+1);
n_cent = n_dir+n_dir/2+1;
gl_trial(n_dir/2+1)=gammalim;
g_trial(n_dir/2+1)=gamma;
for i=1:(n_dir-(2-rem(n_dir,2)))/2
    uip = dir(:,index(n_cent+i));
    angp = acos((ui'*uip)/(norm(ui)*norm(uip)));
    gl_trial(index(n_cent+i)) = (0.5+cos(angp)/2)*gammalim;
    g_trial(index(n_cent+i))  = (0.5+cos(angp)/2)*gamma;
    
    uin = dir(:,index(n_cent-i));
    angn = acos((ui'*uin)/(norm(ui)*norm(uin)));
    gl_trial(index(n_cent-i)) = (0.5+cos(angn)/2)*gammalim;
    g_trial(index(n_cent-i))  = (0.5+cos(angn)/2)*gamma;
end

var.hyp.gammal = max(var.hyp.gammal,gl_trial);
var.hyp.gamma  = max(var.hyp.gamma,g_trial);
var.hyp.gammalh = var.hyp.gammal;
var.hyp.gammah = var.hyp.gamma;

var.est.x=[0;0];
var.est.f=[0;0];
var.est.u=[0;0];
var.est.P=0;
var.est.v=0;
var.hys.uip=0;
var.hys.z=[0;0];
var.hys.phi=0;
var.hys.u=[0;0];

var.G = 1000*1.0238e-04;    
var.A = pi*(var.Ro^2-var.Ri^2)/4;   

var.cav.fc = 3*var.G*var.A;
var.cav.uc = var.cav.fc/var.cav.kv;
var.cav.um = var.cav.uc;
var.cav.fpc=var.cav.fc;
var.cav.upc = var.cav.uc;

var.TSM.th=[0;0];
var.TSM.s=[0;0];
var.TSM.F=[0;0];

var.tol=1e-6;
var.max_iter=100;

%% Run analysis

tux =[0;ux];
tuy =[0;uy];
tP  =[0;P];

F=[0;0];
v=0;
for i=1:length(ux)
    inc=[tux(i+1)-tux(i);tuy(i+1)-tuy(i);tP(i+1)-tP(i)];
    [F(:,i),v(i),var] = HDRB_3D( inc, var );

end

%% Plot figures

figure
plot(sqrt(ux.^2+uy.^2)/var.Ro,-v,sqrt(ux.^2+uy.^2)/var.Ro,-Info(:,3)+7)
ylim([-9.0 0])
grid on


figure
plot(ux/var.h,F(1,:)/A,'r',uy/var.h,F(2,:)/A,'b')
grid on
xlim([-1.8 1.8])
ylim([-0.9 0.9])
xlabel('Strain (\gamma)')
ylabel('Shear stress (MPa)')

