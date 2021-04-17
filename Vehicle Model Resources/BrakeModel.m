function [LinePressure, BreakTorque] = BrakeModel(PedalForce,BalanceBar,BoreDiameter,...
    PadArea,RotorRadius,PedalRatio,PadFriction)

%% Brake Model - Evaluate Brake Forces
% 
% Inputs:
%  PedalForce   - (n,1 numeric) Pedal Force         {F_p} [N-m]
%  BalanceBar   - (n,1 numeric) Percent Balance Bar {%bb} [percent]
%  BoreDiameter - (n,1 numeric) Bore Diameter       {D_b} [m]
%  PadArea      - (n,1 numeric) Pad Area            {A_p} [m^2] 
%  RotorRadius  - (n,1 numeric) Rotor Radius        {} [m]
%  PedalRatio   - (n,1 numeric) Pedal Ratio         {} []
%  PadFriction  - (n,1 numeric) PadFriction         {u_p} 
% 
% Outputs:
%  BreakTorque  - (n,1 numeric) Break Torque        {Tao_b} [N-m]
%  LinePressure - (n,1 numeric) Line Pressure       {P_b} [Pa]
%  
%   
% Notes:
% Test Cases needs to be looked over, commented out parts I'm not sure what 
% Also Not 100% if my units are right in input and output cases
%
% Author(s): 
% Tristan Pham (atlpham@ucdavis.edu) [Sep 2020 - Jun 2021] 

% Last Updated: 16-Apr-2021


%% Test Cases
if nargin == 0
    %%% Test Inputs
    PedalForce = ; 
    
    PercentBalanceBar = ; 
    BoreDiameter = ; 
    
    PadArea = ;
    RotorRadius     = ;
    
    PedalRatio = ;
    PadFriction = ;
    
    fprintf('Executing Brake Model() Test Cases: \n');
    
    
    [LinePressure, BreakTorque] = BrakeModel(PedalForce,PercentBalanceBar,BoreDiameter,...
    PadArea,RotorRadius,PedalRatio,PadFriciton);
    
 % for i = 1:numel()
      %  fprintf('   Steady State Instance %i: \n', i);
      % fprintf('      tau_i = %5.2f [N-m] \n', InputTorque(i));
   % end
    
    
     %return;   
end
    
%% Computation

CylinderForce = [1, 1; (1-BalanceBar), -BalanceBar] \ ...
    [PedalRatio .* PedalForce, 0]';

LinePressure = CylinderForce ./ ( pi/4 .* BoreDiameter.^2 );
BreakTorque = LinePressure/2 .* PadArea .* PadFriction * RotorRadius;