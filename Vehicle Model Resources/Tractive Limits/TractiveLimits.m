clc;clear;close all;
%% Initialize Parameters
%Vehicle
Mass       = 267;                     % Vehicle mass with driver [kg]
TirePos    = [ 0.7625  0.61 0; ...
               0.7625 -0.61 0; ...
              -0.7625  0.61 0; ...
              -0.7625 -0.61 0];       % Contact Patch Location   [m]
Steer      = [0.5; -0.5; -0.5; 0.5];  % Steering Wheel Angle     [degrees]
R_tire     = 0.2032;                  % Radius of tire            [m]
Wheelbase  = 1.525;                   %
TrackWidth = 1.22*ones(1,2);
YawInertia = 130;
CoG        = [(0.5-0.5)*Wheelbase, 0, 0.21];
LatVel     = 0;
YawVel     = 0;

%Tire
EffRadius   = 0.195 .* ones(4,1); %FIXME
Pressure    = 68.9476;
Inclination = [1.5,-1.5,1.5,-1.5]; %FIXME
Idx         = 1;
Model       = struct( 'Pure', 'Pacejka', 'Combined', 'MNC' );
AFx         = -50; %FIXME
AFy         = 0;
AMz         = 0;
load('Hoosier_R25B_16x75-10x7.mat');

%Extra
Downforce = 50; %FIXME
CoP       = [(0.4-0.5)*Wheelbase, 0, 0]; %FIXME
PerLLT    = 0.5; %FIXME 
iter      = 0;
res       = Inf;
tol       = 1e-12;
g         = 9.81; % Gravitational acceleration (m/s^2)


%% Varying Parameter
LongVel = linspace(0,30,100);
for i = 2:length(LongVel)
%% Dependent Parameters
    SpinRate = LongVel(i)/R_tire.* ones(4,1);
    Velocity = LongVel(i);
    NormalLoad_Old = 200 .*ones(4,1);%Initial Normal Load Guess [N]

    %% Slip Estimation
    %Check number of this
    [SlipAngle, SlipRatio, TireVel] = SlipEstimation2( ...
        LongVel(i), LatVel, YawVel, TirePos, Steer, SpinRate, EffRadius);
% SlipRatio = ((SpinRate.*EffRadius)./LongVel(i) - 1);
   while res > tol
    %% Contact Patch Loads
    %Check how many forces we get
        for j = 1:4
            [Fx(j), Fy(j), Mz(j), ~, My(j)] = ContactPatchLoads( Tire, ...
                SlipAngle(j), -SlipRatio(j), ...
                NormalLoad_Old(j), Pressure, Inclination(j), Velocity, ...
                Idx, Model );
        end
    % SlipRatio = linspace(-2,2,100);
    % SlipAngle = 0;
    % for j = 1:length(SlipRatio)
    %     [Fx(j), ~, ~, ~, ~] = ContactPatchLoads( Tire, ...
    %                 SlipAngle, SlipRatio(j), ...
    %                 NormalLoad_Old, Pressure, Inclination, Velocity, ...
    %                 Idx, Model );
    % end
    % Max_Fx = max(Fx);
    
    %         [Fx, Fy, Mz, ~, My] = ContactPatchLoads( Tire, ...
    %             SlipAngle, SlipRatio, ...
    %             NormalLoad_Old, Pressure, Inclination, Velocity, ...
    %             Idx, Model );
    
        %% Vehicle Model
        [LongAcc, LatAcc, YawAcc, LongAccTot, LatAccTot] = ...
            FullTrack3DOFAccelerations( Fx, Fy, Mz, AFx, AFy, AMz, ...
                Wheelbase, TrackWidth, Steer', ...                        
                Mass, YawInertia, CoG, ...                                
                LongVel(i), LatVel, YawVel ) ;
    
        %% Weight Transfer 
        NormalLoad_New = SimplifiedWeightTransfer( LongAccTot, LatAccTot, ...
                Wheelbase, TrackWidth, Mass, CoG, Downforce, CoP, PerLLT ); 
    
        %% Convergence Loop
        res_hold(1) = max(NormalLoad_New(:,:,1) - NormalLoad_Old(1));
        res_hold(2) = max(NormalLoad_New(:,:,2) - NormalLoad_Old(2));
        res_hold(3) = max(NormalLoad_New(:,:,3) - NormalLoad_Old(3));
        res_hold(4) = max(NormalLoad_New(:,:,4) - NormalLoad_Old(4));
        res = max(res_hold);
        NormalLoad_Old(1) = NormalLoad_New(:,:,1);
        NormalLoad_Old(2) = NormalLoad_New(:,:,2);
        NormalLoad_Old(3) = NormalLoad_New(:,:,3);
        NormalLoad_Old(4) = NormalLoad_New(:,:,4);
        iter = iter+1;
    end
    %% Calculate the tractive limit
    mu    = max(Fx)/max(NormalLoad_Old);
    Crr   = max(My)./max(NormalLoad_Old);  %Possibly Fix ME
    h     = 0.21;                                         %CoG Height FIXME
    wr    = NormalLoad_Old(3) + NormalLoad_Old(4);
    %TL(i) = mu*wr*(Wheelbase*0.5-Crr*h)/(Wheelbase-mu*h); %FIXME
    TL(i) = mu*Mass*g*(Wheelbase*0.5-Crr*h)/(Wheelbase+mu*h);
    %TL(i) = (mu*Mass*g*((Wheelbase*0.5-Crr*h)/Wheelbase))/(1-((mu*h)/Wheelbase));
    % TL = (mu * Mass * g ) / (1 + (Crr / mu));
    
    iter = 0;
    res = Inf;
end
%% Generate the graph
plot(round(LongVel(2:end),0), round(TL(2:end),0));
xlabel('Longitudinal Velocity ($\frac{m}{s}$)','Interpreter','latex');
ylabel('Tractive Limit (N)');

%Notes

% TL = (mu * Mass * g * R * (1 - tire_kappa) / (1 + (Crr / mu) ) ) / (1 + (wheel_kappa / (mu / Crr)));
% 
% TL = (mu * Mass * g *R * (1 - tire_kappa) / (1 + (Crr / mu))) / (1 + (wheel_kappa / (mu...
%     / Crr))) - (Tire_Normal_Force * wheel_camber * Sin(delta));
% 
% TL = (mu * Mass * g * R * (1 - tire_kappa)) ./ (1 + (Crr / mu)) - (N * camber * ...
%     sin(steer)) / (1 + (wheel_kappa / (mu / Crr)));
% TL = (mu * Mass * g * R) ./ (1 + (Crr / mu)) - (N * camber * sin(steer)); %N is normal force on single rear tire
% %R is tire radius