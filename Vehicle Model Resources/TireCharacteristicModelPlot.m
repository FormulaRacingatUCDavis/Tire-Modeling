%% Tire Characteristic Model Plotting
clc;clear;close all;

%% Parameters
SlipAngle   = 0;
SlipRatio   = linspace(0,0.2,6);
Velocity    = 10; 
Inclination = 1;
Idx = 1;
NormalLoad  = linspace(200,1000,50);
Pressure    = linspace(62,83,4);
Model       = struct( 'Pure', 'Pacejka', 'Combined', 'MNC' );
%load('Hoosier_R25B_16x75-10x7_Rework.mat');
load('Hoosier_LCO_16x75-10x7.mat');
Pressure_psi = [9,10,11,12];

%% Iterative Variable Initialization
l=0;
k = 0;
q = 1;

%% Plotting For Longitudinal Coeff
figure(1)
for j = 1:4
    for i = 1:6
        
[Fx,Fy,~,~,~] = ContactPatchLoads(Tire,SlipAngle, ...
    SlipRatio(i), NormalLoad, ...
    Pressure(j), Inclination, Velocity, Idx, Model);

%change pressure across row
%change slip ratio across columns 

mux = Fx./NormalLoad;
%muy = Fy./NormalLoad;

l = l+1;
k = k + 1;
subplot(4,6,(l));
plot(NormalLoad,mux,'-');
% plot(NormalLoad,muy,'b','-');
title(sprintf('Pressure = %.1f Psi \n Slip Ratio = %.2f',Pressure_psi(q),SlipRatio(k)))
xlabel('Normal Load(N)')
ylabel('$m_{ux}$','Interpreter','latex')

    end
    k = 0;
    q = q+1;
end
clear;
%% Reinitialize Variables
SlipAngle = linspace(-10,10,6);
SlipRatio = 0;
Velocity = 10; 
Inclination = 1;
Idx = 1;
NormalLoad = linspace(200,1000,50);
Pressure = linspace(62,83,4);
Model       = struct( 'Pure', 'Pacejka', 'Combined', 'MNC' );
%load('Hoosier_R25B_16x75-10x7_Rework.mat');
load('Hoosier_LCO_16x75-10x7.mat');
l=0;
k = 0;
q = 1;
Pressure_psi = [9,10,11,12];

%% Plotting For Lateral Coeff
figure(2)
for j = 1:4
    for i = 1:6
        
[~,Fy,~,~,~] = ContactPatchLoads(Tire,SlipAngle(i), ...
    SlipRatio, NormalLoad, ...
    Pressure(j), Inclination, Velocity, Idx, Model);

%change pressure across row
%change slip ratio across columns 

muy = Fy./NormalLoad;

l = l+1;
k = k + 1;
subplot(4,6,(l));
plot(NormalLoad,muy,'-');
title(sprintf('Pressure = %.1f Psi \n Slip Angle = %.2f',Pressure_psi(q),SlipAngle(k)))
xlabel('Normal Load(N)')
ylabel('$m_{uy}$','Interpreter','latex')

    end
    k = 0;
    q = q+1;
end
