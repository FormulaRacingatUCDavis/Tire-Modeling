function [Variant, Tire] = OverturningVariant(Raw, Tire)
%% Overturning Variant = Fits q and p parameters from Overturning Equation
% Performs Fitting Process on q and p parameters for Overturning Moment
% referenced in Pacejka's "Tire and Vehicle Dynamics" 3rd Edition in
% Section 4.3.2. This will use a contrained fitting process in finding a
% local minima using the fmincon function (although no contraints are used)

% Inputs
%   Raw     - Allocated Data
%   Tire    - Tire Model

% Outputs
%   Variant - Parameters for Overturning Equation
%   Tire    - Updated Tire Model

% Authors
% Carlos Lopez (calopez@ucdavis.edu) [Dec 2020 - June 2022]
% Last Updated: 17- APR - 2021
%% Optimization Variables
qsx1  = optimvar( 'qsx1' , 'Lowerbound',-10   , 'Upperbound', 10    );
qsx2  = optimvar( 'qsx2' , 'Lowerbound',- 5   , 'Upperbound',  5    );
qsx3  = optimvar( 'qsx3' , 'Lowerbound',- 5   , 'Upperbound',  5    );
qsx4  = optimvar( 'qsx4' , 'Lowerbound',- 1   , 'Upperbound',  1    );
qsx5  = optimvar( 'qsx5' , 'Lowerbound',- 5   , 'Upperbound',-0.01  );
qsx6  = optimvar( 'qsx6' , 'Lowerbound',  0   , 'Upperbound',  3    );
qsx7  = optimvar( 'qsx7' , 'Lowerbound',- 5   , 'Upperbound',  5    );
qsx8  = optimvar( 'qsx8' , 'Lowerbound',- 5   , 'Upperbound',  5    );
qsx9  = optimvar( 'qsx9' , 'Lowerbound',- 5   , 'Upperbound',  5    );
qsx10 = optimvar( 'qsx10', 'Lowerbound',  0   , 'Upperbound',  5    );
qsx11 = optimvar( 'qsx11', 'Lowerbound',  0   , 'Upperbound', 10   );

ppmx1 = optimvar( 'ppmx1', 'Lowerbound',-10    , 'Upperbound',10    );

%% Optimization Initialization
x0.qsx1  = -1;
x0.qsx2  = -0.001;
x0.qsx3  = -0.001;
x0.qsx4  = -0.001;
x0.qsx5  = -0.1;
x0.qsx6  = 1;
x0.qsx7  = 0.001;
x0.qsx8  = 0.001;
x0.qsx9  = 0.001;
x0.qsx10 = 0.001;
x0.qsx11 = 0.001;

x0.ppmx1 = 0.001;

%% Optimization Objective 
Obj = fcn2optimexpr(@ErrorMx, qsx1,qsx2,qsx3,qsx4,qsx5,qsx6,qsx7, ...
    qsx8, qsx9, qsx10, qsx11, ppmx1);

%% Optimization Constraint
Constr = optimineq( 0 );

%% Solving Optimization Problem
[Variant.Solution, Variant.Log] = Runfmincon( Obj, x0, Constr, 1);

%% Clearing Optimization Figure
delete( findobj('Type', 'figure', 'Name', 'optimization PlotFcns') );

%% Allocating Solution
Tire.Pacejka.q.S.x(1)  = Variant.Solution.qsx1  ;
Tire.Pacejka.q.S.x(2)  = Variant.Solution.qsx2  ;
Tire.Pacejka.q.S.x(3)  = Variant.Solution.qsx3  ;
Tire.Pacejka.q.S.x(4)  = Variant.Solution.qsx4  ;
Tire.Pacejka.q.S.x(5)  = Variant.Solution.qsx5  ;
Tire.Pacejka.q.S.x(6)  = Variant.Solution.qsx6  ;
Tire.Pacejka.q.S.x(7)  = Variant.Solution.qsx7  ;
Tire.Pacejka.q.S.x(8)  = Variant.Solution.qsx8  ;
Tire.Pacejka.q.S.x(9)  = Variant.Solution.qsx9  ;
Tire.Pacejka.q.S.x(10) = Variant.Solution.qsx10 ;
Tire.Pacejka.q.S.x(11) = Variant.Solution.qsx11 ;

Tire.Pacejka.p.p.mx(1) = Variant.Solution.ppmx1 ;

%% Local Functions
    function [Solution, Log] = Runfmincon( Obj, x0, Constr, n )
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
        
        function stop = OptimLogging( x, optimValues, state )
            stop = false;
            
            if strcmp(state, 'iter')
                Log(ii).x = [Log(ii).x x];
                Log(ii).fval = [Log(ii).fval optimValues.fval];
            end
        end
    end
    
    function RMSE = ErrorMx(qsx1, qsx2, qsx3, qsx4, qsx5, qsx6, qsx7, ...
            qsx8, qsx9, qsx10, qsx11, ppmx1)
        
      [~, Fy, ~, ~, ~] = ContactPatchLoads(Tire, [Raw.Alpha], [Raw.Kappa], ...
          [Raw.Load], [Raw.Pressure], [Raw.Inclination], 10, 1, ...
          struct('Pure', 'Pacejka', 'Combined', 'MNC'));
        
        Mx = Tire.Pacejka.Ro .* [Raw.Load] .* (qsx1 - (qsx2 .* [Raw.Inclination])...
            .* (1 + ppmx1 .* [Raw.dPi]) + qsx3 .* ( Fy./Tire.Pacejka.Fzo)...
            + qsx4 .* cos(qsx5 .* atan(qsx6 .* ([Raw.Load]./Tire.Pacejka.Fzo)).^2) ...
            .* sin(qsx7 .* [Raw.Inclination] + qsx8 .* atan(qsx9 .* (Fy./ ...
            Tire.Pacejka.Fzo))) + qsx10 .* atan(qsx11 .* ([Raw.Load]./Tire.Pacejka.Fzo))...
            .* [Raw.Inclination]);
        RMSE = sqrt( mean( ([Raw.Moment] - Mx).^2) );
    end
end