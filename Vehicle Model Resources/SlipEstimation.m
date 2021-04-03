function [SlipAngle, SlipRatio, TireVel] = SlipEstimation( ...
    LongVel, LatVel, YawVel, TirePos, Steer, SpinRate, EffRadius)

%% SlipEstimation - Slip Angle and Slip Ratio Calculations
% Estimates tire slips and center velocity from chassis velocity and a
% handful of tire states.
% 
% Inputs:
%   LongVel   - (1,1 numeric) Longitudinal Velocity [m/s]
%   LatVel    - (1,1 numeric) Lateral Velocity      [m/s]
%   YawVel    - (1,1 numeric) Yaw Rate              [rad/s]
%   TirePos   - (n,3 numeric) [x,y,0] Tire Position [m]
%   Steer     - (n,1 numeric) Tire Steer Angle      [deg]
%   SpinRate  - (n,1 numeric) Tire Spin Rate        [rad/s]
%   EffRadius - (n,1 numeric) Tire Effective Radius [m]
% 
% Outputs:
%   SlipAngle - (n,1 numeric) Tire Slip Angle      [deg]
%   SlipRatio - (n,1 numeric) Tire Slip Ratio      [ ]
%   TireVel   - (n,1 numeric) Tire Center Velocity [m/s]
%
% Notes:
%   Typically n=4; however, may be vectorized by stacking rows.
%
% Author(s): 
% Blake Christierson (bechristierson@ucdavis.edu) [Sep 2018 - Jun 2021] 
%
% Last Updated: 1-Apr-2021

%%% Test Case
if nargin == 0
    LongVel = 10;
    LatVel = 2;
    YawVel = 0.05;
    TirePos = [0.84 0.66 0; 0.84 -0.66 0; -0.67 0.66 0; -0.67 -0.66 0];
    
    Steer = [25; 20; -0.5; 0.5];
    
    SpinRate = 80 .* ones(4,1);
    EffRadius = 0.19 .* ones(4,1);
    
    warning('Executing SlipEstimation() Test Case')
end

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
