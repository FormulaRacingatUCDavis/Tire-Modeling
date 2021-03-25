function [Variant, Tire] = OverturningVariant(Raw, Tire)
%% Overturning Variant = Fits q and p parameters from Overturning Equation

%% Optimization Variables
qsx1  = optimvar( 'qsx1' , 'Lowerbound',  -15   , 'Upperbound',  15    );
qsx2  = optimvar( 'qsx2' , 'Lowerbound',  -15   , 'Upperbound',  15    );
qsx3  = optimvar( 'qsx3' , 'Lowerbound',  -15   , 'Upperbound',  15    );
qsx4  = optimvar( 'qsx4' , 'Lowerbound',  -15   , 'Upperbound',  15    );
qsx5  = optimvar( 'qsx5' , 'Lowerbound',  -15   , 'Upperbound',  15    );
qsx6  = optimvar( 'qsx6' , 'Lowerbound',  -15   , 'Upperbound',  15    );
qsx7  = optimvar( 'qsx7' , 'Lowerbound',  -15   , 'Upperbound',  15    );
qsx8  = optimvar( 'qsx8' , 'Lowerbound',  -15   , 'Upperbound',  15    );
qsx9  = optimvar( 'qsx9' , 'Lowerbound',  -15   , 'Upperbound',  15    );
qsx10 = optimvar( 'qsx10', 'Lowerbound',  -15   , 'Upperbound',  15    );
qsx11 = optimvar( 'qsx11', 'Lowerbound',  -15   , 'Upperbound',  15    );

ppmx1 = optimvar('ppmx1', 'Lowerbound', -15 , 'Upperbound', 15);

%% Optimization Initialization
x0.qsx1  = 0;
x0.qsx2  = 0;
x0.qsx3  = 0;
x0.qsx4  = 0;
x0.qsx5  = 0;
x0.qsx6  = 0;
x0.qsx7  = 0;
x0.qsx8  = 0;
x0.qsx9  = 0;
x0.qsx10 = 0;
x0.qsx11 = 0;

x0.ppmx1 = 0;

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
        
      [~, Fy, ~, ~, ~] = Tire.ContactPatchLoads([Raw.Alpha], [Raw.Kappa], ...
          [Raw.Load], [Raw.Pressure], [Raw.Inclination], 10, 1, ...
          struct('Pure', 'Pacejka', 'Combined', 'MNC'));
        
        Mx = Tire.Pacejka.Fzo .* [Raw.Load] .* (qsx1 - (qsx2 .* [Raw.Inclination])...
            .* (1 + ppmx1 .* [Raw.dPi]) + qsx3 .* ( Fy./Tire.Pacejka.Fzo)...
            + qsx4 .* cos(qsx5 .* atan(qsx6 .* ([Raw.Load]./Tire.Pacejka.Fzo)).^2) ...
            .* sin(qsx7 .* [Raw.Inclination] + qsx8 .* atan(qsx9 .* (Fy./ ...
            Tire.Pacejka.Fzo))) + qsx10 .* atan(qsx11 .* ([Raw.Load]./Tire.Pacejka.Fzo))...
            .* [Raw.Inclination]);
        RMSE = sqrt( mean( ([Raw.Moment] - Mx).^2) );
    end
end