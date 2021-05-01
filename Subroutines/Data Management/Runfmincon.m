 function [Solution, Log] = Runfmincon( Obj, x0, Constr, n )
%% Overturning Variant = Fits q and p parameters from Overturning Equation
% Performs Fitting Process on q and p parameters for Overturning Moment
% referenced in Pacejka's "Tire and Vehicle Dynamics" 3rd Edition in
% Section 4.3.2. This will use a contrained fitting process in finding a
% local minima using the fmincon function (although no contraints are used)
%
% Inputs:
%   Obj    - Optimization Objective  (1,1 optimexpr)
%   x0     - Initial Solution        (m,1 numeric)
%   Constr - Optimization Constraint (k,1 optimconstr)
%   n      - Number of Multistarts   (1,1 numeric)
%
% Outputs:
%   Solution - Optimization Solution    (m,1 numeric)
%   Log      - Solution & Residual Logs (n,1 struct)
%
% Authors:
% Blake Christierson (bechristierson@ucdavis.edu) [May 2021 - June 2021]
%
% Last Updated: 01-May-2021

%% Problem Setup
% Initializing Logs
Log(n).x = [];
Log(n).fval = [];

% Creating Array of Initial Vectors (Latin Hypercube Sampling)
if nargin > 3 && n > 1
    Jitter = lhsdesign( n-1, numel(fieldnames(x0)) ) - 0.5;

    Initial = (struct2array(x0)/2) .* Jitter + struct2array(x0);
    Initial(:, find( ~struct2array(x0) ) ) = ...
        Jitter( :, find( ~struct2array(x0) ) ) .* 0.2; %#ok<FNDSB>

    Variables = fieldnames( x0 );
    for ii = 1 : n-1
        for jj = 1 : numel( Variables )
            x0(ii+1).(Variables{jj}) = Initial(ii, jj);
        end
    end
end

% Optimization Problem
if isempty(Constr)
    Prob = optimproblem( 'Objective', Obj );
else
    Prob = optimproblem( 'Objective', Obj, 'Constraints', Constr );
end

% Optimization Options
Opts = optimoptions( 'fmincon', ...
    'Algorithm', 'sqp', ...
    'MaxFunctionEvaluations', 10000, ...
    'MaxIterations', 10000, ...
    'Display', 'off', ...
    'PlotFcn', {@optimplotx, @optimplotfval}, ...
    'OutputFcn', @OptimLogging );

% Solving Optimization Problem(s)
for ii = 1 : n
    try
        [Solution(ii), Feval(ii)] = solve( Prob, x0(ii), ...
            'solver', 'fmincon', 'options', Opts ); %#ok<AGROW>
    catch
        continue
    end
end

% Selecting Optimal Solution
[~, MinIdx] = min( Feval );
Solution = Solution(MinIdx);

%% Local Functions 
function stop = OptimLogging( x, optimValues, state )
    stop = false;

    if strcmp(state, 'iter')
        Log(ii).x = [Log(ii).x x];
        Log(ii).fval = [Log(ii).fval optimValues.fval];
    end
end

end
