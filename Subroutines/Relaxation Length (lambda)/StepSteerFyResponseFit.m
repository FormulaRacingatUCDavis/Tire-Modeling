function [FitResult, GoF, RL] = StepSteerFyResponseFit( t, Fy, V)
%% Relaxation Length Curve Fitting - Transient Tire Model
% This fits filtered raw tire data to obtain tau, which is used to
% calculate relaxation length.
%  
% Inputs:
%   Tire           - (struct)      Pacejka Parameters
%   SlipAngle      - (n,1 numeric) Slip Angle          {alpha} [deg]
%   NormalLoad     - (n,1 numeric) Normal Load         {F_z}   [N]
%   LateralForce   - (n,2 numeric) Lateral Force       {F_y}   [N]
%   AligningMoment - (n,3 numeric) Aligning Moment     {M_z}   [Nm]
%   Pressure       - (n,1 numeric) Inflation Pressure  {P_i}   [kPa]
%   Velocity       - (n,1 numeric) Center Velocity     {v_c}   [m/s]
%   Model          - (struct)      Fidelity Choices 
%
% Outputs:
%   RelaxationLength - (n,1 numeric) Relaxation Length {lambda} [m]
%
% Notes:
%   2nd order curve fitting is currently not implemented (6/2/21) because
%   accuracy of relaxation length data is currently sufficient to justify 
%   tire choice.
%
% Author(s): 
% Blake Christierson (bechristierson@ucdavis.edu) [Sep 2018 - Jun 2021] 
% Leonardo Howard    (leohoward@ucdavis.edu     ) [Feb 2021 -         ]
% 
% Last Updated: 02-June-2021

%% Initialization.

% Initialize arrays to store fits and goodness-of-fit.
FitResult = cell( 2, 1 );
GoF = struct( 'sse', cell( 2, 1 ), ...
    'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );

%% Fit: 'First Order'.
[xData, yData] = prepareCurveData( t, Fy );

% Set up fittype and options.
ft = fittype( 'A*(1-exp(-t/tau))', 'independent', 't', 'dependent', 'Fy' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf 0];
opts.StartPoint = [1 0.05];
opts.Upper = [Inf Inf];

% Fit model to data.
[FitResult{1}, GoF(1)] = fit( xData, yData, ft, opts );

%% First Order RL calculation.
CurveFitValues = coeffvalues(FitResult{1}); %extract coefficients from the curve fitting
Tau = CurveFitValues(2); %extract tau from the curve fitting coefficients
        
RL = Tau.*mode(V); %obtain the relaxation length values  

%% Fit: 'Second Order'.
%{
[xData, yData] = prepareCurveData( t, Fy );

% Set up fittype and options.
ft = fittype( 'K*(1 - sqrt((zeta*w0)^2 + w0*sqrt(1-zeta^2))/(w0*sqrt(1-zeta^2))*exp(-zeta*w0*t)*sin(w0*sqrt(1-zeta^2)*t+acos(zeta)))', 'independent', 't', 'dependent', 'Fy' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf 0 0];
opts.StartPoint = [1 10 0.5];
opts.Upper = [Inf Inf 1];

% Fit model to data.
[FitResult{2}, GoF(2)] = fit( xData, yData, ft, opts );
%}
