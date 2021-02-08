clc; clear; close all;

%% Friction Ellipse Plotting
% Plots friction ellipse of a desired tire

%% Load Tire Model
load('Hoosier_R25B_16x75-10x7.mat')

%% Operating Condition
Pi  = 12 ;
Fz  = 675;
Inc = 0  ;
Vx  = 10 ;

Fidelity.Pure     = 'Pacejka';
Fidelity.Combined = 'MNC'    ;

%% Generate Slip Angle Isocontours
Alpha = linspace(-15,15,20);
Kappa = linspace(-0.15, 0.15, 100);

[Alpha, Kappa] = meshgrid( Alpha, Kappa );

[Fx, Fy, ~, ~, ~] = Tire.ContactPatchLoads( Alpha, Kappa, Fz, Pi, Inc, Vx, 1, Fidelity );

plot( Fy, Fx, 'k' ); hold on;

%% Generate Slip Angle Isocontours
Kappa = linspace(-0.15,0.15,20);
Alpha = linspace(-15, 15, 100);

[Kappa, Alpha] = meshgrid( Kappa, Alpha );

[Fx, Fy, ~, ~, ~] = Tire.ContactPatchLoads( Alpha, Kappa, Fz, Pi, Inc, Vx, 1, Fidelity );

plot( Fy, Fx, 'k' ); axis equal;

ylim( [-1.1.*max(abs(Fx), [], 'all') 1.1.*max(abs(Fx), [], 'all')] )
xlim( [-1.1.*max(abs(Fy), [], 'all') 1.1.*max(abs(Fy), [], 'all')] )

ylabel( 'Longitudinal Force: $F_{x}$ [$N$]' )
xlabel( 'Lateral Force: $F_{y}$ [$N$]' )