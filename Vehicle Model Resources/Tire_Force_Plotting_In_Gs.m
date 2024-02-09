clc;clear;close all;

Mass = 264.5/4; %Mass of (Car+driver)/4
g = 9.81;
load('Hoosier_R25B_16x75-10x7.mat');
Pressure    = 70;
Inclination = 1;
Velocity    = 10;
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
    
    Fx = Fx/(Mass*g);
    Fy = Fy/(Mass*g);

    figure(1)
    sgtitle( {'Slip-Load Surfaces', ...
        ['$P_{i} ='  , num2str(Pressure)   , '$ [$kPa$],' ...
          '$\gamma = ', num2str(Inclination), '$ [$deg$]']} ,'Interpreter','latex')
      
    subplot(1,2,1)
    surf( SlipRatio, NormalLoad, Fx )
    xlabel( '$\kappa$ [ ]' ,'Interpreter','latex')
    ylabel( '$F_{z}$ [$N$]','Interpreter','latex' )
    zlabel( '$F_{x}$ [$g$]' ,'Interpreter','latex')
    title( 'Longitudinal Force' )
    
    subplot(1,2,2)
    surf( SlipAngle, NormalLoad, Fy )
    xlabel( '$\alpha$ [$deg$]','Interpreter','latex' )
    ylabel( '$F_{z}$ [$N$]','Interpreter','latex' )
    zlabel( '$F_{y}$ [$g$]','Interpreter','latex' )
    title( 'Lateral Force' )

    figure(2)
    subplot(1,2,1)
    hold on
    for i = 1:25
        plot(NormalLoad(i,:),Fy(i,:))
    end
    title("Tire Lateral Force vs Normal Load w/ Increasing Slip Angles")
    xlabel("Normal Load (N)")
    ylabel("Lateral Force (g)")
    
    subplot(1,2,2)
    hold on
    for i = 1:51
        dy = diff(Fy(i,:))./diff(NormalLoad(i,:));
        plot(NormalLoad(i,2:end),dy)
    end
    title("Lateral Force Sensitivity vs Normal Load w/ Increasing Slip Angles")
    xlabel("Normal Load (N)")
    ylabel("Lateral Force Sensitivity")

%% Regular Plotting
clear;
Mass = 264.5/4; %Mass of (Car+driver)/4
g = 9.81;
load('Hoosier_R25B_16x75-10x7.mat');
Pressure    = linspace(62,83,4);
Inclination = 1.7;
Velocity    = 10;
Idx         = 1;
Model       = struct( 'Pure', 'Pacejka', 'Combined', 'MNC' );
SlipRatio = linspace(- 1, 1,51);
%SlipAngle = linspace(-20,20,51);
SlipAngle = 0;
NormalLoad  = linspace(200,1000,5);
Pressure_psi = [9,10,11,12];
% Iterative Variable Initialization
l=0;
k = 0;
q = 1;

% Plotting For Longitudinal
figure(3)
sgtitle( 'Longitudinal Forces vs Slip Ratio','Interpreter','latex')
for j = 1:4
    for i = 1:5

        for p = 1:length(SlipRatio)
[Fx(p),~,~,~,~] = ContactPatchLoads(Tire,SlipAngle, ...
    SlipRatio(p), NormalLoad(i), ...
    Pressure(j), Inclination, Velocity, Idx, Model);
        end

%change pressure across row
%change slip ratio across columns 

Fx = (3/4).*Fx./(Mass*g);
%muy = Fy./NormalLoad;

l = l+1;
k = k + 1;
subplot(4,5,(l));
plot(SlipRatio,Fx,'-');
% plot(NormalLoad,muy,'b','-');
title(sprintf('Pressure = %.1f Psi \n Normal Load = %.2f',Pressure_psi(q),NormalLoad(k)))
xlabel('SlipRatio ($deg$)','Interpreter','latex')
ylabel('$Fx (g)$','Interpreter','latex')

    end
    k = 0;
    q = q+1;
end

%% Regular Plotting
clear;
Mass = 264.5/4; %Mass of (Car+driver)/4
g = 9.81;
load('Hoosier_R25B_16x75-10x7.mat');
Pressure    = linspace(62,83,4);
Inclination = 0;
Velocity    = 10;
Idx         = 1;
Model       = struct( 'Pure', 'Pacejka', 'Combined', 'MNC' );
SlipRatio = 0;
SlipAngle = linspace(-20,20,51);
NormalLoad  = linspace(200,1000,5);
Pressure_psi = [9,10,11,12];

% Iterative Variable Initialization
l=0;
k = 0;
q = 1;

% Plotting For Longitudinal
figure(4)
sgtitle( 'Lateral Forces vs Slip Angle','Interpreter','latex')
for j = 1:4
    for i = 1:5

        for p = 1:length(SlipAngle)
            [~,Fy(p),~,~,~] = ContactPatchLoads(Tire,SlipAngle(p), ...
                SlipRatio, NormalLoad(i), ...
                Pressure(j), Inclination, Velocity, Idx, Model);
        end

%change pressure across row
%change slip ratio across columns 

Fy = Fy./(Mass*g);
%muy = Fy./NormalLoad;

l = l+1;
k = k + 1;
subplot(4,5,(l));
plot(SlipAngle,Fy,'-');
% plot(NormalLoad,muy,'b','-');
title(sprintf('Pressure = %.1f Kpa \n Normal Load = %.2f',Pressure_psi(q),NormalLoad(k)))
xlabel('SlipAngle ($deg$)','Interpreter','latex')
ylabel('$Fy (g)$','Interpreter','latex')

    end
    k = 0;
    q = q+1;
end

%%
clear;
Mass = 301.99/4; %Mass of (Car+driver)/4
g = 9.81;
load('Hoosier_R25B_16x75-10x7.mat');
Pressure    = 70;
Inclination = 1;
Velocity    = 10;
Idx         = 1;
Model       = struct( 'Pure', 'Pacejka', 'Combined', 'MNC' );

%%% Slip-Load Plots 
%     SlipRatio = linspace(- 1, 1,51);
%     SlipAngle = linspace(-20,20,51);
SlipAngle = -2;
SlipRatio = 0;
    
    NormalLoad  = linspace(0,800,20);
        
%     [~         , SlipRatio] = meshgrid( NormalLoad, SlipRatio );
%     [NormalLoad, SlipAngle] = meshgrid( NormalLoad, SlipAngle );
    
    [Fx, Fy, Mz, Mx, My] = ContactPatchLoads( Tire, ...
        SlipAngle, SlipRatio, ...
        NormalLoad, Pressure, Inclination, Velocity, ...
        Idx, Model );
    
    Fx = Fx/(Mass*g);
    Fy = Fy/(Mass*g);

Fy_Combined = fliplr(Fy) + Fy;
NormalLoad = NormalLoad * 2;
figure(5)
plot(NormalLoad,Fy_Combined)
title("Pairs Grip Analysis")
xlabel("Normal Load For Single Tire(N)")
ylabel("Combined Grip of Front Axel (g)")