function [ Response ] = PureLateralResponseSurfaces( Raw, Mesh, Nominal, Tire )

%% Defining Operating Condition Functions
dPi = @(Pi) (Pi - Tire.Pacejka.Pio) ./ Tire.Pacejka.Pio;
dFz = @(Fz) (Fz - Tire.Pacejka.Fzo) ./ Tire.Pacejka.Fzo;

%% Fitting C_{y} Surface
% Initial Vector
x0.pcy1 = mean( [Nominal.C] );

C.Surface = @(Pi, Fz, Gam, x0) x0.pcy1 + 0*(Pi + Fz + Gam);

%% Fitting D_{y} Surface
% Initial Vector 
D0.pdy1 = 2;
D0.pdy2 = -0.5;
D0.pdy3 = 0;
D0.ppy3 = 0;
D0.ppy4 = 0;

% Optimization Variables
pdy1 = optimvar( 'pdy1', 'Lowerbound', 0 );
pdy2 = optimvar( 'pdy2', 'Upperbound', 0 );
pdy3 = optimvar( 'pdy3' );

ppy3 = optimvar( 'ppy3' );
ppy4 = optimvar( 'ppy4' );

% Optimization Objective
Obj = fcn2optimexpr( @ErrorDy, pdy1, pdy2, pdy3, ppy3, ppy4 );

% Solving Optimization Problem
[D.Solution, D.Log] = Runfmincon( Obj, D0, [], 3 );
 
% Allocating Optimal Solution
x0.pdy1 = D.Solution.pdy1;
x0.pdy2 = D.Solution.pdy2;
x0.pdy3 = D.Solution.pdy3;

x0.ppy3 = D.Solution.ppy3;
x0.ppy4 = D.Solution.ppy4;

D.Surface = @(Pi, Fz, Gam, x0) ( x0.pdy1 + x0.pdy2.*dFz(Fz) ) .* ...
            ( 1 + x0.ppy3.*dPi(Pi) + x0.ppy4.*dPi(Pi).^2 ) .* ...
            ( 1 - x0.pdy3.*Gam.^2 ) .* Fz;

clear pdy1 pdy2 pdy3 ppy3 ppy4 D0

%% Fitting E_{y} Surface
% Initial Vector 
E0.pey1 = 0;
E0.pey2 = 0;
E0.pey5 = -0.01;

% Optimization Variables
pey1 = optimvar( 'pey1', 'UpperBound', 1 );
pey2 = optimvar( 'pey2' );
pey5 = optimvar( 'pey5', 'UpperBound', 0 );

% Optimization Objective
Obj = fcn2optimexpr( @ErrorEy, pey1, pey2, pey5 );

% Optimization Constraint
Constr = pey1 - pey2 <= 1;

% Solving Optimization Problem
[E.Solution, E.Log] = Runfmincon( Obj, E0, Constr, 3 );

% Allocating Solution
x0.pey1 = E.Solution.pey1;
x0.pey2 = E.Solution.pey2;
x0.pey3 = 0;
x0.pey4 = 0;
x0.pey5 = E.Solution.pey5;

E.Surface = @(Pi, Fz, Gam, x0) (x0.pey1 + x0.pey2.*dFz(Fz)) .* ...
    (1 + x0.pey5.*Gam.^2);

clear pey1 pey2 pey5 E0

%% Fitting K_{y \alpha} Surface
% Initial Vector 
K0.pky1 = -4;
K0.pky2 = 0.25;
K0.pky3 = 0;
K0.pky4 = 2;
K0.pky5 = 0;

K0.ppy1 = 0;
K0.ppy2 = 0;

% Optimization Variables
pky1 = optimvar( 'pky1', 'Upperbound',-0.1 );
pky2 = optimvar( 'pky2', 'Lowerbound', 0.1, 'Upperbound', 5   );
pky3 = optimvar( 'pky3', 'Lowerbound',-5  , 'Upperbound', 5   );
pky4 = optimvar( 'pky4', 'Lowerbound', 1.5, 'Upperbound', 2.5 );
pky5 = optimvar( 'pky5', 'Lowerbound',-1  , 'Upperbound', 1   );

ppy1 = optimvar( 'ppy1', 'Lowerbound',-1, 'Upperbound', 1 );
ppy2 = optimvar( 'ppy2', 'Lowerbound',-1, 'Upperbound', 1 );

% Optimization Objective
Obj = fcn2optimexpr( @ErrorKya, pky1, pky2, pky3, pky4, pky5, ppy1, ppy2 );

% Solving Optimization Problem
[K.Solution, K.Log] = Runfmincon( Obj, K0, [], 10 );

% Allocating Solution
x0.pky1 = K.Solution.pky1;
x0.pky2 = K.Solution.pky2;
x0.pky3 = K.Solution.pky3;
x0.pky4 = K.Solution.pky4;
x0.pky5 = K.Solution.pky5;

x0.ppy1 = K.Solution.ppy1;
x0.ppy2 = K.Solution.ppy2;

K.Surface = @(Pi, Fz, Gam, x0) ( x0.pky1 .* Tire.Pacejka.Fzo .* ( 1 + x0.ppy1.*dPi(Pi) ) .* ...
    ( 1 - x0.pky3.*abs(Gam) ) .* sin( x0.pky4 .* atan( (Fz./Tire.Pacejka.Fzo) ./ ...
    ( ( x0.pky2 + x0.pky5.*Gam.^2 ) .* ( 1 + x0.ppy2.*dPi(Pi) ) ) ) ) ) ;

clear pky1 pky2 pky3 pky4 pky5 ppy1 ppy2 K0

%% Fitting V_{y} Surface
% Initial Vector 
V0.pvy1 = 0;
V0.pvy2 = 0;
V0.pvy3 = 0;
V0.pvy4 = 0;

% Optimization Variables
pvy1 = optimvar( 'pvy1' );
pvy2 = optimvar( 'pvy2' );
pvy3 = optimvar( 'pvy3' );
pvy4 = optimvar( 'pvy4' );

% Optimization Objective
Obj = fcn2optimexpr( @ErrorVy, pvy1, pvy2, pvy3, pvy4 );

% Solving Optimization Problem
[V.Solution, V.Log] = Runfmincon( Obj, V0, [], 3 );

% Allocating Solution
x0.pvy1 = V.Solution.pvy1;
x0.pvy2 = V.Solution.pvy2;
x0.pvy3 = V.Solution.pvy3;
x0.pvy4 = V.Solution.pvy4;

V.Surface = @(Pi, Fz, Gam, x0) Fz .* ( ( x0.pvy1 + x0.pvy2.*dFz(Fz) ) + ...
    ( x0.pvy3 + x0.pvy4.*dFz(Fz) ).*Gam );

clear pvy1 pvy2 pvy3 pvy4 V0

%% Fitting K_{yg0} Surface
% Initial Vector
Kyg0.pky6 = 0;
Kyg0.pky7 = 0;

Kyg0.ppy5 = 0;

% Optimization Variables
pky6 = optimvar( 'pky6', 'Lowerbound', -5, 'Upperbound', 5 );
pky7 = optimvar( 'pky7', 'Lowerbound', -5, 'Upperbound', 5 );

ppy5 = optimvar( 'ppy5', 'Lowerbound', -5, 'Upperbound', 5 );

% Optimization Objective
Obj = fcn2optimexpr( @ErrorKyg, pky6, pky7, ppy5 );

% Solving Optimization Problem
[Kyg.Solution, Kyg.Log] = Runfmincon( Obj, Kyg0, [], 3 ); 

% Allocating Solution
x0.pky6 = Kyg.Solution.pky6;
x0.pky7 = Kyg.Solution.pky7;

x0.ppy5 = Kyg.Solution.ppy5;

%% Fitting H_{y} Surface
% Initial Vector 
H0.phy1 = 0;
H0.phy2 = 0;

% Optimization Variables
phy1 = optimvar( 'phy1', 'Lowerbound', -5, 'Upperbound', 5 );
phy2 = optimvar( 'phy2', 'Lowerbound', -5, 'Upperbound', 5 );

% Optimization Objective
Obj = fcn2optimexpr( @ErrorHy, phy1, phy2 );

% Solving Optimization Problem
[H.Solution, H.Log] = Runfmincon( Obj, H0, [], 3 ); 

% Allocating Solution
x0.phy1 = H.Solution.phy1;
x0.phy2 = H.Solution.phy2;

H.Surface = @(Pi, Fz, Gam, x0) ( x0.phy1 + x0.phy2.*dFz(Fz) ) + ...
    Fz.*Gam .* ( ( x0.pky6 + x0.pky7.*dFz(Fz) ) .* ( 1 + x0.ppy5.*dPi(Pi) ) - ...
    ( x0.pvy3 + x0.pvy4.*dFz(Fz) ) ) ./ Surface.K(Pi, Fz, Gam);

clear phy1 phy2 pky6 pky7 ppy5 H0

%% Alpha-Sort of x0 
x0 = orderfields( x0, sort(fieldnames( x0 )) );

%% Clearing Optimization Figure
delete( findobj( 'Type', 'figure', 'Name', 'Optimization PlotFcns' ) );

%% Allocating Outputs
Response.x0 = x0;

Response.C = C;
Response.D = D;
Response.E = E;
Response.K = K;
Response.H = H;
Response.V = V;

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

    function MeanSquareError = ErrorDy( pdy1, pdy2, pdy3, ppy3, ppy4 )
        DySurface = (pdy1 + pdy2.*[Mesh.dFz]) .* ...
            (1 + ppy3.*[Mesh.dPi] + ppy4.*[Mesh.dPi].^2) .* ...
            (1 - pdy3.*[Mesh.Inclination].^2).*[Mesh.Load];
        
        MeanSquareError = mean( ( [Nominal.D] - DySurface ).^2 );
    end

    function MeanAbsoluteError = ErrorEy( pey1, pey2, pey5 )
        EySurface = (pey1 + pey2.*[Mesh.dFz]) .* ...
            (1 + pey5.*[Mesh.Inclination].^2);
        
        MeanAbsoluteError = mean( abs( [Nominal.E] - EySurface ) );
    end

    function MeanSquareError = ErrorKya( pky1, pky2, pky3, pky4, pky5, ppy1, ppy2 )
        KySurface = pky1 .* Tire.Pacejka.Fzo .* ( 1 + ppy1.*[Mesh.dPi] ) .* ...
            ( 1 - pky3.*abs([Mesh.Inclination]) ) .* sin( pky4 .* ...
            atan( ([Mesh.Load]./Tire.Pacejka.Fzo) ./ ...
            ( ( pky2 + pky5.*[Mesh.Inclination].^2 ) .* ( 1 + ppy2.*[Mesh.dPi] ) ) ) );
        
        MeanSquareError = mean( ( [Nominal.K] - KySurface ).^2 );
    end

    function MeanSquareError = ErrorVy( pvy1, pvy2, pvy3, pvy4 )
        VySurface = [Mesh.Load] .* ( ( pvy1 + pvy2.*[Mesh.dFz] ) + ...
            ( pvy3 + pvy4.*[Mesh.dFz] ).*[Mesh.Inclination] );
        
        MeanSquareError = mean( ( [Nominal.V] - VySurface ).^2 );
    end
    
    function MeanSquareError = ErrorKyg( pky6, pky7, ppy5 )
        Data.Slip        = [Raw.Slip];
        Data.Load        = [Raw.Load];
        Data.dFz         = [Raw.dFz];
        Data.Inclination = [Raw.Inclination];
        Data.dPi         = [Raw.dPi];
        Data.Force       = [Raw.Force];
        
        Data.Load        = Data.Load(        abs(Data.Slip) < deg2rad(0.2) );
        Data.dFz         = Data.dFz(         abs(Data.Slip) < deg2rad(0.2) );
        Data.Inclination = Data.Inclination( abs(Data.Slip) < deg2rad(0.2) );
        Data.dPi         = Data.dPi(         abs(Data.Slip) < deg2rad(0.2) );
        Data.Force       = Data.Force(       abs(Data.Slip) < deg2rad(0.2) );
        
        KygSurface = Data.Load.*Data.Inclination .* ...
            ( pky6 + pky7.*Data.dFz ) .* ( 1 + ppy5.*Data.dPi );
        
        MeanSquareError = mean( ( Data.Force - KygSurface ).^2 );
    end

    function MeanSquareError = ErrorHy( phy1, phy2 )
        HySurface = ( phy1 + phy2.*[Mesh.dFz] ) + [Mesh.Load].*[Mesh.Inclination] .* ...
            ( ( x0.pky6 + x0.pky7.*[Mesh.dFz] ) .* ( 1 + x0.ppy5.*[Mesh.dPi] ) - ...
            ( x0.pvy3 + x0.pvy4.*[Mesh.dFz] ) ) ./ ...
            K.Surface( [Mesh.Pressure], [Mesh.Load], [Mesh.Inclination], x0 );
        
        MeanSquareError = mean( ( [Nominal.H] - HySurface ).^2 );
    end
end