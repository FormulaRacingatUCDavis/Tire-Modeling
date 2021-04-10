function [InputTorque, SpinAcc] = WheelSpeed( SpinRate, ...
    DriveTorque, BrakeTorque, TractiveForce, EffRadius, Inertia, Damping)

%% WheelSpeed - Wheel Torque Balance & Spin Acceleration
% Calculates brake or drive torque required for steady state applications
% or the spin acceleration for transient applications.
% 
% Inputs:
%   SpinRate      - (n,1 numeric) Wheel Spin Rate       {omega} [rad/s]
%   DriveTorque   - (n,1 numeric) Drive Torque          {tau_D} [N-m]
%   BrakeTorque   - (n,1 numeric) Brake Torque          {tau_B} [N-m]
%   TractiveForce - (n,1 numeric) Tractive Force        {F_x}   [N]
%   EffRadius     - (n,1 numeric) Tire Effective Radius {r_e}   [m]
%   Inertia       - (n,1 numeric) Spin Inertia          {I_s}   [kg-m^2]
%   Damping       - (n,1 numeric) Spin Damping          {b_s}   [N-m-s/rad]
% 
% Outputs:
%   InputTorque - (n,1 numeric) Axle Input Torque      {tau_{D|B}} [N-m]
%   SpinAcc     - (n,1 numeric) Tire Spin Acceleration {omega_dot} [rad/s^2]
%
% Notes:
%   Steady state computation will be enacted if BrakeTorque or DriveTorque
%   arguments are left empty.
%
% Author(s): 
% Blake Christierson (bechristierson@ucdavis.edu) [Sep 2018 - Jun 2021] 
%
% Last Updated: 7-Apr-2021

%%% Test Cases
if nargin == 0
    %%% Test Inputs
    SpinRate = 80;
    
    DriveTorque = 20;
    BrakeTorque = 0;
    
    TractiveForce = 300;
    EffRadius     = 0.19;
    
    Inertia = 
    Damping = 
end

end

