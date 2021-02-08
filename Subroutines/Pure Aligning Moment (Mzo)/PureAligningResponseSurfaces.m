function [ Response ] = PureAligningResponseSurfaces( Mesh, Nominal, Tire )

%% Defining Operating Condition Functions
dPi = @(Pi) (Pi - Tire.Pacejka.Pio) ./ Tire.Pacejka.Pio;
dFz = @(Fz) (Fz - Tire.Pacejka.Fzo) ./ Tire.Pacejka.Fzo;

%% Fitting B_{t} Surface
% Preparing Data
[Bt.Data(:,1), Bt.Data(:,2), Bt.Data(:,3)] = ...
    prepareSurfaceData( [Mesh.dFz], [Mesh.Inclination], [Nominal.Bt] );

% Set up fittype and options.
Bt.Type = fittype( ['(qbz1 + qbz2*dFz + qbz3*dFz^2) .*' ...
    '(1 + qbz5*abs(Gam) + qbz6*Gam^2)'], ...
    'independent', {'dFz', 'Gam'}, 'dependent', 'Bt' );

Bt.Opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
Bt.Opts.Robust = 'LAR';

Bt.Opts.Display = 'Off';

Bt.Opts.MaxFunEvals = 10000;
Bt.Opts.MaxIter = 10000;

Bt.Opts.StartPoint = [0 0 0 0 0];
Bt.Opts.Lower = [-5 -5 -5 -5 -5];
Bt.Opts.Upper = [5 5 5 5 5];

% Fit model to data.
[Bt.Solution, Bt.Gof] = fit( [Bt.Data(:,1), Bt.Data(:,2)], Bt.Data(:,3), Bt.Type, Bt.Opts );
 
% Allocating Optimal Solution
x0.qbz1 = Bt.Solution.qbz1;
x0.qbz2 = Bt.Solution.qbz2;
x0.qbz3 = Bt.Solution.qbz3;

x0.qbz5 = Bt.Solution.qbz5;
x0.qbz6 = Bt.Solution.qbz6;

%% Fitting C_{t} Surface
% Initial Vector
x0.qcz1 = mean( [Nominal.Ct] );

Ct.Surface = @(Pi, Fz, Gam, x0) x0.qcz1 + 0*(Pi + Fz + Gam);

%% Fitting D_{t} Surface
% Initial Vector 
Dt0.qdz1 = 1;
Dt0.qdz2 = -1;
Dt0.qdz3 = 0;
Dt0.qdz4 = 0;

Dt0.ppz1 = 0;

% Optimization Variables
qdz1 = optimvar( 'qdz1' );
qdz2 = optimvar( 'qdz2' );
qdz3 = optimvar( 'qdz3' );
qdz4 = optimvar( 'qdz4' );

ppz1 = optimvar( 'ppz1' );

% Optimization Objective
Obj = fcn2optimexpr( @ErrorDt, qdz1, qdz2, qdz3, qdz4, ppz1 );

% Solving Optimization Problem
[Dt.Solution, Dt.Log] = Runfmincon( Obj, Dt0, [], 3 );
 
% Allocating Optimal Solution
x0.qdz1 = Dt.Solution.qdz1;
x0.qdz2 = Dt.Solution.qdz2;
x0.qdz3 = Dt.Solution.qdz3;
x0.qdz4 = Dt.Solution.qdz4;

x0.ppz1 = Dt.Solution.ppz1;
        
Dt.Surface = @(Pi, Fz, Gam, x0) Tire.Pacejka.Ro .* (Fz./Tire.Pacejka.Fzo) .* ...
    ( x0.qdz1 + x0.qdz2.*dFz(Fz) ) .* ( 1 - x0.ppz1.*dPi(Pi) ) .* ...
    ( 1 + x0.qdz3.*abs(Gam) + x0.qdz4.*Gam.^2 );


clear qdz1 qdz2 qdz3 qdz4 ppz1 Dt0

%% Fitting H_{t} 
% Preparing Data
[Ht.Data(:,1), Ht.Data(:,2), Ht.Data(:,3)] = ...
    prepareSurfaceData( [Mesh.dFz], [Mesh.Inclination], [Nominal.Ht] );

% Set up fittype and options.
Ht.Type = fittype( 'qhz1 + qhz2*dFz + (qhz3 + qhz4*dFz)*Gam', ...
    'independent', {'dFz', 'Gam'}, 'dependent', 'Ht' );

Ht.Opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
Ht.Opts.Robust = 'LAR';

Ht.Opts.Display = 'Off';

Ht.Opts.MaxFunEvals = 10000;
Ht.Opts.MaxIter = 10000;

Ht.Opts.StartPoint = [0 0 0 0];
Ht.Opts.Lower = [-5 -5 -5 -5];
Ht.Opts.Upper = [5 5 5 5];

% Fit model to data.
[Ht.Solution, Ht.Gof] = fit( [Ht.Data(:,1), Ht.Data(:,2)], Ht.Data(:,3), Ht.Type, Ht.Opts );
 
% Allocating Optimal Solution
x0.qhz1 = Ht.Solution.qhz1;
x0.qhz2 = Ht.Solution.qhz2;
x0.qhz3 = Ht.Solution.qhz3;
x0.qhz4 = Ht.Solution.qhz4;

clear qhz1 qhz2 qhz3 qhz4 Ht0

%% Fitting E_{t}
% Initial Vector 
Et0.qez1 = 0.5;
Et0.qez2 = 0;
Et0.qez3 = 0;
Et0.qez4 = 0;
Et0.qez5 = 0;

% Optimization Variables
qez1 = optimvar( 'qez1', 'Lowerbound', -5, 'Upperbound', 5 );
qez2 = optimvar( 'qez2', 'Lowerbound', -5, 'Upperbound', 5 );
qez3 = optimvar( 'qez3', 'Lowerbound', -5, 'Upperbound',-0.001);
qez4 = optimvar( 'qez4', 'Lowerbound', -5, 'Upperbound', 5 );
qez5 = optimvar( 'qez5', 'UpperBound', 0 );

% Optimization Objective
Obj = fcn2optimexpr( @ErrorEt, qez1, qez2, qez3, qez4, qez5 );

% Solving Optimization Problem
[Et.Solution, Et.Log] = Runfmincon( Obj, Et0, [], 3 );

% Allocating Solution
x0.qez1 = Et.Solution.qez1;
x0.qez2 = Et.Solution.qez2;
x0.qez3 = Et.Solution.qez3;
x0.qez4 = Et.Solution.qez4;
x0.qez5 = Et.Solution.qez5;

Et.Surface = @(Pi, Fz, Gam, x0) (x0.qez1 + x0.qez2.*dFz(Fz) + x0.qez3.*dFz(Fz).^2) .* ...
    (1 + (x0.qez4 + x0.qez5.*Gam).*(2/pi).* ...
    atan( Bt.Solution( dFz(Fz), Gam ).*x0.qcz1.*5 ) );

clear qez1 qez2 qez3 qez4 qez5 Et0

%% Fitting B_{r} Surface
% Initial Vector
x0.qbz10 = mean( [Nominal.qbz10] );

%% Fitting D_{r} Surface
% Initial Vector 
Dr0.qdz6  = 0;
Dr0.qdz7  = 0;
Dr0.qdz8  = 0;
Dr0.qdz9  = 0;
Dr0.qdz10 = 0;
Dr0.qdz11 = 0;

Dr0.ppz2  = 0;

% Optimization Variables
qdz6 = optimvar( 'qdz6' );
qdz7 = optimvar( 'qdz7' );
qdz8 = optimvar( 'qdz8' );
qdz9 = optimvar( 'qdz9' );
qdz10 = optimvar( 'qdz10' );
qdz11 = optimvar( 'qdz11' );

ppz2 = optimvar( 'ppz2' );

% Optimization Objective
Obj = fcn2optimexpr( @ErrorDr, qdz6, qdz7, qdz8, qdz9, qdz10, qdz11, ppz2 );

% Solving Optimization Problem
[Dr.Solution, Dr.Log] = Runfmincon( Obj, Dr0, [], 3 );
 
% Allocating Optimal Solution
x0.qdz6 = Dr.Solution.qdz6;
x0.qdz7 = Dr.Solution.qdz7;
x0.qdz8 = Dr.Solution.qdz8;
x0.qdz9 = Dr.Solution.qdz9;
x0.qdz10 = Dr.Solution.qdz10;
x0.qdz11 = Dr.Solution.qdz11;

x0.ppz2 = Dr.Solution.ppz2;
        
Dr.Surface = @(Pi, Fz, Gam, x0)  Tire.Pacejka.Ro .* Fz .* ( ( x0.qdz6 + x0.qdz7.*dFz(Fz) ) + ...
            ( ( x0.qdz8 + x0.qdz9.*dFz(Fz) ) .* ( 1 + x0.ppz2.*dPi(Pi) ) + ...
            ( x0.qdz10 + x0.qdz11.*dFz(Fz) ) .* abs(Gam) ).*Gam );

clear qdz6 qdz7 qdz8 qdz9 qdz10 qdz11 ppz2 Dr0

%% Alpha-Sort of x0 
x0 = orderfields( x0, sort(fieldnames( x0 )) );

%% Clearing Optimization Figure
delete( findobj( 'Type', 'figure', 'Name', 'Optimization PlotFcns' ) );

%% Allocating Outputs
Response.x0 = x0;

Response.Bt = Bt;
Response.Ct = Ct;
Response.Dt = Dt;
Response.Et = Et;
Response.Ht = Ht;

%Response.Br = Br;
Response.Dr = Dr;

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

    function RMSE = ErrorDt( qdz1, qdz2, qdz3, qdz4, ppz1 )
        DtSurface = Tire.Pacejka.Ro .* ([Mesh.Load]./Tire.Pacejka.Fzo) .* ...
            ( qdz1 + qdz2.*[Mesh.dFz] ) .* ( 1 - ppz1.*[Mesh.dPi] ) .* ...
            ( 1 + qdz3.*abs([Mesh.Inclination]) + qdz4.*[Mesh.Inclination].^2 );
        
        RMSE = sqrt( mean( ([Nominal.Dt] - DtSurface).^2 ) );
    end

    function RMSE = ErrorEt( qez1, qez2, qez3, qez4, qez5 )
        EtSurface = (qez1 + qez2.*[Mesh.dFz] + qez3.*[Mesh.dFz].^2) .* ...
            (1 + (qez4 + qez5.*[Mesh.Inclination]).*(2/pi).* ...
            atan( Bt.Solution([Mesh.dFz], [Mesh.Inclination]).*x0.qcz1.*5 ) );
        
        RMSE = sqrt( mean( ([Nominal.Et] - EtSurface).^2 ) );
    end
    
    function RMSE = ErrorDr( qdz6, qdz7, qdz8, qdz9, qdz10, qdz11, ppz2 )
        DrSurface = Tire.Pacejka.Ro .* [Mesh.Load] .* ( ( qdz6 + qdz7.*[Mesh.dFz] ) + ...
            ( ( qdz8 + qdz9.*[Mesh.dFz] ) .* ( 1 + ppz2.*[Mesh.dPi] ) + ...
            ( qdz10 + qdz11.*[Mesh.dFz] ) .* abs([Mesh.Inclination]) ).*[Mesh.Inclination] );
        
        RMSE = sqrt( mean( ([Nominal.Dr] - DrSurface).^2 ) );
    end
end