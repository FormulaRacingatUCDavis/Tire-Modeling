%aligning moment (Mz), SlipAngle, NormalLoad plot (3D)
clear;close all;clc;
%load('Hoosier_R25B_16x75-10x7_Rework.mat');
load('Hoosier_LCO_16x75-10x7.mat');
SlipAngle = linspace(-20,20,50);
SlipRatio = linspace(-1,1,50);
Velocity = 10; 
Inclination = 1;
Idx = 1;
NormalLoad = linspace(0,2000,50);
Pressure =  70; %linspace(69,90,1);
Model       = struct( 'Pure', 'Pacejka', 'Combined', 'MNC' );

[X,Y] = meshgrid(SlipAngle,NormalLoad);
[~,~,Mz,~,~] = ContactPatchLoads(Tire,X,SlipRatio,Y, ...
    Pressure,Inclination,Velocity,Idx,Model);
surf(X,Y,Mz)

title('Aligning Moment vs Slip Angle vs Normal Load')
a=num2str(Pressure);
b=num2str(Velocity);
c=num2str(Inclination);
e='At ';
g='kPa, ';
h='m/s, ';
i=' deg of Inclination';
subtitle(append(e,a,g,b,h,c,i));
xlabel('Slip Angle');
ylabel('Normal Load');
zlabel('Aligning Moment (Mz)');