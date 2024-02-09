%camberthrust
clear;clc;close all;
%load('Hoosier_R25B_16x75-10x7_Rework.mat');
load('Hoosier_LCO_16x75-10x7.mat');

points = 50;
SlipAngle = linspace(-10,10,points);
SlipRatio = linspace(-1,1,points);
Velocity = 0; 
Inclination = linspace(-4,4,points);
Idx = 1;
NormalLoad = linspace(200, 2000, points);%linspace(0,2000,50);
Pressure =  70; 
Model       = struct( 'Pure', 'Pacejka', 'Combined', 'MNC' );

%%surf1 = inclination vs normal load vs Fy[lateralforce])
 [X1,Y1] = meshgrid(Inclination, NormalLoad);

 [Fx,Fy,~,~,~] = ContactPatchLoads(Tire,SlipAngle,SlipRatio,Y1, ...
     Pressure,X1,Velocity,Idx,Model);
 figure(1);
 surf(X1,Y1,Fy)
xlabel('Inclination');
ylabel('Normal Load');
zlabel('Lateral Force');
 
%%surf2 = mux vs Fx vs NormalLoad
%mux=Fx/NormalLoad;
%Inclination = 1;
%mux = reshape(mux',1,[]);
%[X2,Y2] = meshgrid(mux,NormalLoad);
%[Fx,Fy,~,~,~] = ContactPatchLoads(Tire,SlipAngle,SlipRatiFo,Y, ...
 %   Pressure,Inclination,Velocity,Idx,Model);
%figure(2);
%surf(X2,Y2,Fy)
%xlabel('mux');
%ylabel('Normal Load');
%zlabel('Fy');
%surf2 = muy vs Fy vs NormalLoad
%muy=Fy./NormalLoad;
%muy = reshape(muy',1,[]);
%[X3,Y3] = meshgrid(muy, NormalLoad);
%[Fx,Fy,~,~,~] = ContactPatchLoads(Tire,SlipAngle,SlipRatio,Y, ...
%    Pressure,X,Velocity,Idx,Model);
%figure(3);
%surf(X3,Y3,Fy)
%xlabel('muy');
%ylabel('Normal Load');
%zlabel('Fy');