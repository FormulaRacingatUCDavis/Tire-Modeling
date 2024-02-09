clear;close all;clc;
%load('Hoosier_R25B_16x75-10x7_Rework.mat');
load('Hoosier_LCO_16x75-10x7.mat');
SlipAngle = linspace(-20,20,50);
SlipRatio = linspace(-1,1,50);
Velocity = 0; %linspace(0,40,50); 
Inclination = 1;
Idx = 1;
NormalLoad = linspace(200, 2000, 12);%linspace(0,2000,50);
Pressure =  70; %linspace(69,90,1);
Model       = struct( 'Pure', 'Pacejka', 'Combined', 'MNC' );

%[X,Y] = meshgrid(SlipAngle,NormalLoad);
Fx=[];
Fy=[];
for k=1:12
for i= 1:50
    for j= 1:50
[Fx(i,j),Fy(i,j),~,~,~] = ContactPatchLoads(Tire,SlipAngle(i),SlipRatio(j),NormalLoad(k), ...
    Pressure,Inclination,Velocity,Idx,Model);
    end
end

subplot(4, 3, k)
plot(Fx, Fy)
a=num2str(NormalLoad(k));
b=('Normal Load: ');
title(append(b,a))
xlabel('$Lateral Force$', 'Interpreter', 'latex')
ylabel('$Longitudinal Force$','Interpreter', 'latex')
end