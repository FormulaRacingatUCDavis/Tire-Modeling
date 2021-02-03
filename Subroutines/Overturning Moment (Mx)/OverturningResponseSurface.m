function [ Response ] = OverturningResponseSurface(Mesh, Nominal, Tire)

%% Defining Operating Condition Functions

dPi = @(Pi) (Pi - Tire.Pacejka.Pio)./ Tire.Pacejka.Pio;
dFz = @(Fz) (Fz - Tire.Pacejka.Fzo)./ Tire.Pacejka.Fzo;

%% Fitting Mx Surface
% Initial Vector
Mx0.qsx1  = 0;
Mx0.qsx2  = 0;
Mx0.qsx3  = 0;
Mx0.qsx4  = 0;
Mx0.qsx5  = 0;
Mx0.qsx6  = 0;
Mx0.qsx7  = 0;
Mx0.qsx8  = 0;
Mx0.qsx9  = 0;
Mx0.qsx10 = 0;
Mx0.qsx11 = 0;

Mx0.ppMx1 = 0;

% Optimization Variables
qsx1  = optimvar( 'qsx1'  );
qsx2  = optimvar( 'qsx2'  );
qsx3  = optimvar( 'qsx3'  );
qsx4  = optimvar( 'qsx4'  );
qsx5  = optimvar( 'qsx5'  );
qsx6  = optimvar( 'qsx6'  );
qsx7  = optimvar( 'qsx7'  );
qsx8  = optimvar( 'qsx8'  );
qsx9  = optimvar( 'qsx9'  );
qsx10 = optimvar( 'qsx10' );
qsx11 = optimvar( 'qsx11' );

ppMx1 = optimvar( 'ppMx1' );

% Optimization Objective
Obj = fcn2optimexpr( @ErrorMx, qsx1, qsx2, qsx3, qsx4, qsx5, qsx6, qsx7, ...
    qsx8, qsx9, qsx10, qsx11, ppMx1);

% Solving Optimization Problem
[Mx.Solution, Mx.Log] = Runfmincon(Obj, Mx, [], 3);

% Allocating Solution
x0.qsx1  = Mx.Solution.qsx1  ;
x0.qsx2  = Mx.Solution.qsx2  ;
x0.qsx3  = Mx.Solution.qsx3  ;
x0.qsx4  = Mx.Solution.qsx4  ;
x0.qsx5  = Mx.Solution.qsx5  ;
x0.qsx6  = Mx.Solution.qsx6  ;
x0.qsx7  = Mx.Solution.qsx7  ;
x0.qsx8  = Mx.Solution.qsx8  ;
x0.qsx9  = Mx.Solution.qsx9  ;
x0.qsx10 = Mx.Solution.qsx10 ;
x0.qsx11 = Mx.Solution.qsx11 ;

x0.ppMx1 = Mx.Solution.ppMx1 ;


Mx.Surface = @(Fz, Gam, Pi, x0) (Tire.Pacejka.Ro * Tire.Pacejka.Fzo) * ... 
    ( x0.qsx1 - ((x0.qsx2 * Gam) * (1 + x0.ppMx1 * dPi(Pi))) +(x0.qsx3 * ...
    (Fy(Fz, Gam, Pi, x0)./Tire.Pacejka.Fzo)) + (x0.qsx4 * cos(x0.qsx5 * ...
    atan(x0.qsx6 * (Fz./Tire.Pacejka.Fzo).^2))* sin((x0.qsx7 * Gam) + ...
    (x0.qsx8 * atan( x0.qsx9 * (Fy(Fz,Gam,Pi,x0)./Tire.Pacejka.Fzo))))) + ... 
    (x0.qsx10 * atan(x0.qsx11 * (Fz./Tire.Pacejka.Fzo)) * Gam));

clear qsx1 qsx2 qsx3 qsx4 qsx5 qsx6 qsx7 qsx8 qsx9 qsx10 qsx11 ppMx1

%% Alpha-Sort of x0
x0 = orderfields( x0, sort(fieldnames( x0 )) );

%% Clearing Optimization Figure
delete( findobj( 'Type', 'figure', 'Name', 'Optimization PlotFcns' ) );

%% Allocating Outputs
Response.x0 = x0;

Response.Mx = Mx;

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
    function MeanSquareError = ErrorMx(qsx1, qsx2, qsx3, qsx4, qsx5, ...
            qsx6, qsx7, qsx8, qsx9, qsx10, qsx11, ppMx1)
        MxSurface = (Tire.Pacejka.Ro * Tire.Pacejka.Fzo) * ( qsx1 - ... 
            ((qsx2 * [Mesh.Inclination]) * (1 + ppMx1 * [Mesh.dPi])) +(qsx3 * ...
            (Fy./Tire.Pacejka.Fzo)) + (qsx4 * cos(qsx5 * ...
            atan(qsx6 * ([Mesh.Load]./Tire.Pacejka.Fzo).^2))* sin((qsx7 * [Mesh.Inclination]) + ...
            (qsx8 * atan( qsx9 * (Fy./Tire.Pacejka.Fzo))))) + ... 
            (qsx10 * atan(qsx11 * ([Mesh.Load]./Tire.Pacejka.Fzo)) * [Mesh.Inclination]));
        MeanSquareError = mean( ([Raw.Moment]  - MxSurface).^2);
    end
end
