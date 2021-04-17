function [SlipAngle, SlipRatio, TireVel] = SlipEstimation( ...
    LongVel, LatVel, YawVel, TirePos, Steer, SpinRate, EffRadius)

%% SlipEstimation - Slip Angle and Slip Ratio Calculations
% Estimates tire slips and center velocity from chassis velocity and a
% handful of tire states.
% 
% Inputs:
%   LongVel   - (1,1 numeric) Longitudinal Velocity {x_dot}   [m/s]
%   LatVel    - (1,1 numeric) Lateral Velocity      {y_dot}   [m/s]
%   YawVel    - (1,1 numeric) Yaw Rate              {psi_dot} [rad/s]
%   TirePos   - (n,3 numeric) [x,y,0] Tire Position {p}       [m]
%   Steer     - (n,1 numeric) Tire Steer Angle      {delta}   [deg]
%   SpinRate  - (n,1 numeric) Tire Spin Rate        {omega}   [rad/s]
%   EffRadius - (n,1 numeric) Tire Effective Radius {r_e}     [m]
% 
% Outputs:
%   SlipAngle - (n,1 numeric) Tire Slip Angle      {alpha} [deg]
%   SlipRatio - (n,1 numeric) Tire Slip Ratio      {kappa} [ ]
%   TireVel   - (n,1 numeric) Tire Center Velocity {v_c}   [m/s]
%
% Notes:
%   Typically n=4; however, may be vectorized by stacking rows.
%
% Author(s): 
% Blake Christierson (bechristierson@ucdavis.edu) [Sep 2018 - Jun 2021] 
%
% Last Updated: 1-Apr-2021

%% Test Cases
if nargin == 0
    LongVel = 10;
    LatVel = 2;
    YawVel = 0.1;
    TirePos = [ 0.84  0.66 0; ...
                0.84 -0.66 0; ...
               -0.67  0.66 0; ...
               -0.67 -0.66 0];
    
    Steer = [25; 20; -0.5; 0.5];
    
    SpinRate = 100 .* ones(4,1);
    EffRadius = 0.19 .* ones(4,1);
    
    [SlipAngle, SlipRatio, TireVel] = SlipEstimation( ...
        LongVel, LatVel, YawVel, TirePos, Steer, SpinRate, EffRadius);

    fprintf('Executing SlipEstimation() Test Cases: \n');
    for i = 1:numel(SlipAngle)
        fprintf('   Instance %i: \n', i);
        fprintf('      alpha = %5.2f [deg] \n', SlipAngle(i));
        fprintf('      kappa = %5.2f [ ] \n', SlipRatio(i));
        fprintf('      v_c   = %5.2f [m/s] \n', TireVel(i));
    end
    
    return;
end

%% Computation
%%% Number of Evaluations
n = numel(Steer);

%%% Tire Center Velocities
TireVel = zeros(n,3);
for i = 1:n
    TireVel(i,:) = [LongVel; LatVel; 0] + cross( [0; 0; YawVel], TirePos(i,:)' );
end

%%% Slip Angle
SlipAngle = atan2d( TireVel(:,2), TireVel(:,1) ) - Steer;

%%% Slip Ratio
SlipRatio = ( SpinRate .* EffRadius ) ./ ( norm(TireVel,2) .* cosd( SlipAngle ) ) - 1;
