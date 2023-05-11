clear 
close all

%% Load of experimental results and model parameters

load('../Data/Test.mat')
param=importdata('../Data/Parameters.txt');

var.hyp.a1= param(1);
var.hyp.a2= param(2);
var.hyp.a3= param(3);
var.hyp.fs1= param(4);
var.hyp.ps1= param(5);
var.hyp.fs2= param(6);
var.hyp.ps2= param(7);
var.hyp.fs3= param(8);
var.hyp.ps3= param(9);
var.hyp.fm = param(10);
var.hyp.pm = param(11);
var.hys.uy= param(12);
var.hys.fy= param(13);
var.hys.B= param(14);
var.hys.n= param(15);
var.hys.phimax= param(16);
var.hys.pphi= param(17);

var.hyp.h = h;

%% Internal variables

% Modify this in case of the device has been tested previously and the
% initial variables has changed

var.hys.k = var.hys.fy/var.hys.uy;
var.hys.gamma = 1-var.hys.B;

var.hyp.gp=0;
var.hyp.gn=0;
var.hyp.gcp=0;
var.hyp.gcn=0;
var.hyp.glp=0;
var.hyp.gln=0;
var.hys.ui=0;
var.hys.uip=0;
var.hys.z=0;
var.hys.phi=0;
var.hys.up=0;
var.est.x=0;
var.est.f=0;

%% Solution

fm=[0;0;0];xm=0;km=0;
for j=2:length(d)
    inc=d(j)-d(j-1);
    [fm(:,j),xm(j),km(j),var] = URB(inc,var);
end

%% Figure

fig = figure;
plot(xm,fm(3,:),'r',d,f,'k',LineWidth=1.5)
title('Comparisson between numerical and experimental results')
xlabel('Displacement [cm]')
ylabel('Force [tonf]')
grid on
legend('Proposed','Experimental','Location','NorthWest')
saveas(fig,'../Figures/Matlab.png')