% clc; clear; close all
format short g

%%
warning('Executing SlipEstimation() Test Case')

addpath( genpath( fileparts( which( 'ContactPatchLoads.m' ) ) ) );
%load('Models\TestTire.mat'); %#ok<LOAD> missing Tire.Pacejka.p.P.Mx
%load('Hoosier_R25B_18x75-10x7.mat');

%%% Old Tire: Hoosier_R25B_16x75-10x7.mat
% load('Hoosier_R25B_16x75-10x7.mat');
%load('TESTING-1-19-2025-R25B-Default.mat')
%load('Tire.mat');

% %%% New Tires as of 2024: Hoosier_43075_16x75-10_R20-8.mat
% load('Hoosier_43075_16x75-10_R20-8.mat')
 load('Hoosier_R20_16(18)x75(60)-10x8(7).mat')

%%% Nominal Test Case Conditions
Pressure    = 70;
Inclination = 1;
%Velocity    = 10;
Velocity    = 0;
Idx         = 1;
Model       = struct( 'Pure', 'Pacejka', 'Combined', 'MNC' );

%%% Slip-Load Plots 
SlipRatio = linspace(- 1, 1,51);
SlipAngle = linspace(-20,20,51);

NormalLoad  = 0:100:2000;
    
[~         , SlipRatio] = meshgrid( NormalLoad, SlipRatio );
[NormalLoad, SlipAngle] = meshgrid( NormalLoad, SlipAngle );

[Fx, Fy, Mz, Mx, My] = ContactPatchLoads( Tire, ...
    SlipAngle, SlipRatio, ...
    NormalLoad, Pressure, Inclination, Velocity, ...
    Idx, Model );
 
figure
sgtitle( {'Slip-Load Surfaces', ...
    ['$P_{i} ='  , num2str(Pressure)   , '$ [$kPa$],' ...
      '$\gamma = ', num2str(Inclination), '$ [$deg$]']} ,'Interpreter','latex')
  
subplot(3,2,1)
surf( SlipRatio, NormalLoad, Fx )
xlabel( '$\kappa$ [ ]' ,'Interpreter','latex')
ylabel( '$F_{z}$ [$N$]','Interpreter','latex' )
zlabel( '$F_{x}$ [$N$]' ,'Interpreter','latex')
title( 'Longitudinal Force' )

subplot(3,2,2)
surf( SlipAngle, NormalLoad, Fy )
xlabel( '$\alpha$ [$deg$]','Interpreter','latex' )
ylabel( '$F_{z}$ [$N$]','Interpreter','latex' )
zlabel( '$F_{y}$ [$N$]','Interpreter','latex' )
title( 'Lateral Force' )

subplot(3,2,3)
surf( SlipAngle, NormalLoad, Mz )
xlabel( '$\alpha$ [$deg$]' ,'Interpreter','latex')
ylabel( '$F_{z}$ [$N$]','Interpreter','latex' )
zlabel( '$M_{z}$ [$Nm$]','Interpreter','latex' )
title( 'Aligning Moment' )

subplot(3,2,4)
surf( SlipAngle, NormalLoad, Mx )
xlabel( '$\alpha$ [$deg$]' ,'Interpreter','latex')
ylabel( '$F_{z}$ [$N$]' ,'Interpreter','latex')
zlabel( '$M_{x}$ [$Nm$]','Interpreter','latex' )
title( 'Overturning Moment' )

subplot(3,2,5)
surf( SlipRatio, NormalLoad, My )
xlabel( '$\kappa$ [ ]','Interpreter','latex' )
ylabel( '$F_{z}$ [$N$]','Interpreter','latex' )
zlabel( '$M_{y}$ [$Nm$]','Interpreter','latex' )
title( 'Rolling Resistance' )

% subplot(3,2,6)
% Fz = linspace(0,2000,51)';
% plot(Fy/Fz, SlipAngle)
  
%%% Friction Ellipse Plotting
SlipRatio = linspace(-1, 1,51);
SlipAngle = linspace(-20,20,51);

[SlipRatio, SlipAngle] = meshgrid( SlipRatio, SlipAngle );

NormalLoad = 700;

[Fx, Fy, ~, ~, ~] = ContactPatchLoads( Tire, ...
    SlipAngle, SlipRatio, ...
    NormalLoad, Pressure, Inclination, Velocity, ...
    Idx, Model );

ColorMap = colormap( 'parula' );
figure
for i = 1 : size( SlipAngle, 1)
    plot( Fy(i,:), Fx(i,:), 'Color', ...
        ColorMap( round( (SlipAngle(i,1)-SlipAngle(1)) .* ...
            (length(ColorMap)-1) ./ diff( SlipAngle([1 end]) )) + 1, : ) ); hold on;
end

for j = 1 : size( SlipRatio, 2)
    plot( Fy(:,j), Fx(:,j), 'Color', ...
        ColorMap( round( (SlipRatio(1,j)-SlipRatio(1)) .* ...
            (length(ColorMap)-1) ./ diff( SlipRatio([1 end]) )) + 1, : ) );
end
xlabel( 'Lateral Force, $F_{y}$ [$N$]','Interpreter','latex' ); 
ylabel( 'Longitudinal Force, $F_{x}$ [$N$]','Interpreter','latex' );
title( 'Friction Ellipse' )
ylim([-2000, 2000]);
xlim([-2000, 2000]);
grid on


%%% Inclination Sensitivity
Velocity = 15;
Inclination = -3:0.05:3;
NormalLoad = 700; %Newtons

SlipRatio = linspace(- 1, 1,51);
SlipAngle = linspace(-20,20,51);

[~          , SlipRatio] = meshgrid( Inclination, SlipRatio );
[Inclination, SlipAngle] = meshgrid( Inclination, SlipAngle );

[Fx, Fy, ~, ~, ~] = ContactPatchLoads( Tire, ...
    SlipAngle, SlipRatio, ...
    NormalLoad, Pressure, Inclination, Velocity, ...
    Idx, Model );

figure
hold on;
plot( Inclination(1,:), max(abs(Fx),[],1)./max(abs(Fx),[],'all') );
plot( Inclination(1,:), max(abs(Fy),[],1)./max(abs(Fy),[],'all') )
xlabel( 'Inclination, $\gamma$ [$deg$]','Interpreter','latex' )
ylabel( 'Normalized Grip' )
title( 'Inclination Sensitivity')
legend( '$F_{x}$', '$F_{y}$','Interpreter','latex' )

figure
hold on;
plot( Inclination(1,:), max(abs(Fx),[],1)./NormalLoad );
plot( Inclination(1,:), max(abs(Fy),[],1)./NormalLoad )
xlabel( 'Inclination, $\gamma$ [$deg$]','Interpreter','latex' )
ylabel( 'Normalized Grip' )
title( 'Inclination Sensitivity')
legend( '$F_{x}$', '$F_{y}$','Interpreter','latex' )

figure
hold on;
surf(Inclination, SlipAngle, Fy./NormalLoad)
grid on
view(3)
xlabel( 'Inclination, $\gamma$ [$deg$]','Interpreter','latex' )
ylabel( 'Slip Angle (deg)' ,'Interpreter','latex' )
zlabel( 'Lateral Force Friction Coefficient','Interpreter','latex' )



% 
% %%% Coefficient of Friction
% Pressure    = 70;
% Inclination = 1;
% % Velocity    = 10;
% Velocity = 0;
% Idx         = 1;
% Model       = struct( 'Pure', 'Pacejka', 'Combined', 'MNC' );
% 
% SlipRatio = linspace(0, 1,51);
% SlipAngle = linspace(0,20,51);
% NormalLoad  = 0:100:2000;
% 
% [~         , SlipRatio] = meshgrid( NormalLoad, SlipRatio );
% [NormalLoad, SlipAngle] = meshgrid( NormalLoad, SlipAngle );
% 
% [Fx, Fy, ~,~,~] = ContactPatchLoads( Tire, ...
%     SlipAngle, SlipRatio, ...
%     NormalLoad, Pressure, Inclination, Velocity, ...
%     Idx, Model );
% 
% muY = -Fy ./ NormalLoad;
% 
% figure 
% hold on
% grid()
% surf(SlipAngle, NormalLoad, muY);
% view(3)
% xlabel("Slip Angle [deg]")
% ylabel("Normal Load [N]")
% zlabel("\mu_y [-]")
% title("Coefficient of Friction in Y")
% 
% Mat.NormalLoad = NormalLoad(find(NormalLoad == 700));
% Mat.SlipAngle = SlipAngle(find(NormalLoad == 700));
% Mat.muY = muY(find(NormalLoad == 700));
% 
% figure
% hold on
% grid()
% plot(Mat.SlipAngle, Mat.muY)
% xlabel("Slip Angle [deg]")
% ylabel("\mu_y [-]")
% title("Coefficient of Friction in Y at 700 N")
% 
% 
% muX = Fx ./ NormalLoad;
% 
% figure 
% hold on
% grid()
% surf(SlipRatio, NormalLoad, muX);
% view(3)
% xlabel("Slip Ratio ([0-1]",'Interpreter','latex');
% ylabel("Normal Load $(F_z)$ [N]",'Interpreter','latex');
% zlabel("Coefficient of Friction $(\mu_x)$ [-]",'Interpreter','latex');
% title("Coefficient of Friction in X",'Interpreter','latex');
% 
% Mat.NormalLoad = NormalLoad(find(NormalLoad == 700));
% Mat.SlipRatio = SlipRatio(find(NormalLoad == 700));
% Mat.muX = muX(find(NormalLoad == 700));
% 
% figure
% hold on
% grid()
% plot(Mat.SlipRatio, Mat.muX)
% xlabel("Slip Ratio [0-1]")
% ylabel("\mu_x [-]")
% title("Coefficient of Friction in X at 700 N")
% 





%%
% %%% Single Evaluation
% Inclination = 1;
% SlipAngle = 5;
% SlipRatio = 0.03;
% 
% [Fx, Fy, Mz, Mx, My] = ContactPatchLoads( Tire, ...
%     SlipAngle, SlipRatio, ...
%     NormalLoad, Pressure, Inclination, Velocity, ...
%     Idx, Model );