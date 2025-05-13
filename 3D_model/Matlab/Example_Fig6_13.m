clear
clc
close all

%% Input

th=180*pi/180;

%% Load patern
n_dir=72;
x=[0:0.01:1 0.99:-0.01:0.01];
x=[x x x x x x];
np=length(x);
x0=zeros(1,np);

d1=[x0 0.25*x 0.5*x x0 x0 0.25*x 0.5*x x x0 x0 x0 0.25*x 0.5*x x 1.5*x x0 x0 x0 x0 0.25*x 0.5*x x 1.5*x 2.5*x x0 x0 x0 x0 x0];
d2=[0.25*x x0 x0 0.25*x 0.5*x x0 x0 x0 0.25*x 0.5*x x x0 x0 x0 x0 0.25*x 0.5*x x 1.5*x x0 x0 x0 x0 x0 0.25*x 0.5*x x 1.5*x 2.5*x];

list_changes = np*[1 3 5 8 11 15 19 24];
a=find(d1~=0);
b=find(d2~=0);
d=[d1+cos(th)*d2;d2*sin(th)];

%% Parameters

var.hys.k=8.0;
var.hys.fy=0.2;
var.hys.uy=var.hys.fy/var.hys.k;
var.hys.B=0.6;
var.hys.gamma=1-var.hys.B;
var.hys.n=0.312;
var.hys.phimax=0.183; 
var.hys.pphi=0.008;

var.hyp.a1= 1.911607967;
var.hyp.a2= 0.17486988;
var.hyp.a3= 0.280000;
var.hyp.fs1= 0.6400000; 
var.hyp.ps1= 1.3000000;
var.hyp.fs2= 15.457; 
var.hyp.ps2= 0.304; 
var.hyp.fs3= 0.5000000;
var.hyp.ps3= 0.235; 
var.hyp.fm = 0.165; 
var.hyp.pm = 0.7000000;
var.hyp.fc = 0.8;
var.hyp.h  =1;
var.h  = 1;

%% Update variables

angle=(0:360/n_dir:360*(n_dir-1)/n_dir)*pi/180;
var.hyp.dir=[cos(angle);sin(angle)];
var.hyp.gamma  =zeros(1,n_dir);
var.hyp.gammal =zeros(1,n_dir);
var.hyp.gammac =zeros(1,n_dir);
var.hyp.gammalh = var.hyp.gammal;
var.hyp.gammah = var.hyp.gamma;

%% Run Analysis

var.est.x=[0;0];
var.est.f=[0;0];
var.est.u=[0;0];
var.est.P=0;
var.est.v=0;
var.hys.uip=0;
var.hys.z=[0;0];
var.hys.phi=0;
var.hys.u=[0;0];

for i = 2:length(d)
    
    inc=d(:,i)-d(:,i-1);
    [fs(:,i),var] = BIHDRB3( inc,var );
   
    % Comment for figure 13
    if rem(i,np)==0
        var.est.u=[0;0];
        var.hyp.gammac =zeros(1,n_dir);
        var.hys.z=[0;0];
    end


end

%% Plot Results

figr = figure;
plot(d1(a),fs(1,a))
hold on
plot(d2(b),fs(1,b)*cos(-th)-fs(2,b)*sin(-th), '--')
xlim([-0.25 2.75])
ylim([-0.5 4.5])
grid on
% exportgraphics(figr,'figures/Result_'+string(th*180/pi)+'.pdf','ContentType','vector')