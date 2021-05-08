function [InputTorque, SpinAcc] = WheelSpeed( SpinRate,DriveTorque, ...
    BrakeTorque, TractiveForce, EffRadius, Inertia, Damping )
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
% Joseph Sanchez     (jomsanchez@ucdavis.edu)     [Sep 2020 - Jun 2022]

% Last Updated: 10-Apr-2021

%% Test Cases
if nargin == 0
    %%% Test Inputs
    SpinRate = 20:20:100; 
    
    DriveTorque = 20*ones(size(SpinRate)); 
    BrakeTorque = 0*ones(size(SpinRate)); 
    
    TractiveForce = 300*ones(size(SpinRate));
    EffRadius     = 0.19*ones(size(SpinRate));
    
    Inertia = 0.148;
    Damping = 0.01;
    
    fprintf('Executing WheelSpeed() Test Cases: \n');
    
    %%% Steady State Test Cases
    [InputTorque, ~] = WheelSpeed( SpinRate, [], [], ...
        TractiveForce, EffRadius, Inertia, Damping );
    
    for i = 1:numel(InputTorque)
        fprintf('   Steady State Instance %i: \n', i);
        fprintf('      tau_i = %5.2f [N-m] \n', InputTorque(i));
    end
    
    %%% Transient Test Cases
    [~, SpinAcc] = WheelSpeed( SpinRate, DriveTorque, BrakeTorque, ...
        TractiveForce, EffRadius, Inertia, Damping );  
    
    for i = 1:numel(SpinAcc)
        fprintf('   Transient Instance %i: \n', i);
        fprintf('      omega_dot = %5.2f [rad/s^2] \n', SpinAcc(i));
    end
    
     return;   
end

%% Computation
if isempty( DriveTorque ) || isempty( BrakeTorque )
    %%% Steady State Input Torque
    InputTorque = TractiveForce.*EffRadius + Damping.*SpinRate;
        
    SpinAcc = NaN;
else
    %%% Transient Spin Acceleration
    SpinAcc = ( DriveTorque - TractiveForce.*EffRadius - BrakeTorque - ...
        Damping.*SpinRate ) ./ Inertia;
    
    InputTorque = NaN;
end

