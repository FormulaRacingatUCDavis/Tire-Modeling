function [ Response ] = PureLongitudinalResponseSurfaces( Tire, Mesh, Nominal )
%% PureLongitudinalResponseSurfaces - Fits Parameter Responses to Nominal Results
% Fits response surfaces for C,D,K,E,H,V to initialize full variant pure
% slip curve fit.
% 
% Inputs:
%   Tire    - Tire Model
%   Mesh    - Experimental Operating Conditions
%   Nominal - Fitted P6 Parameters
%
% Inputs:
%   Response - Fitted Response Surface Parameters
%
% Author(s): 
% Blake Christierson (bechristierson@ucdavis.edu) [Sep 2018 - Jun 2021] 
% Carlos Lopez       (calopez@ucdavis.edu       ) [Jan 2019 -         ]
% 
% Last Updated: 15-Feb-2021

%% Defining Operating Condition Functions
dPi = @(Pi) (Pi - Tire.Pacejka.Pio) ./ Tire.Pacejka.Pio;
dFz = @(Fz) (Fz - Tire.Pacejka.Fzo) ./ Tire.Pacejka.Fzo;

%% Fitting C_{x} Surface
% Initial Vector
x0.pcx1 = mean( [Nominal.C] );

C.Surface = @(Pi, Fz, Gam, x0) x0.pcx1 + 0*(Pi + Fz + Gam);

%% Fitting D_{x} Surface
% Initial Vector 
D0.pdx1 = 2;
D0.pdx2 = -0.5;
D0.pdx3 = 0;
D0.ppx3 = 0;
D0.ppx4 = 0;

% Optimization Variables
pdx1 = optimvar( 'pdx1', 'Lowerbound', 0 );
pdx2 = optimvar( 'pdx2', 'Upperbound', 0 );
pdx3 = optimvar( 'pdx3' );

ppx3 = optimvar( 'ppx3' );
ppx4 = optimvar( 'ppx4' );

% Optimization Objective
Obj = fcn2optimexpr( @ErrorDx, pdx1, pdx2, pdx3, ppx3, ppx4 );

% Solving Optimization Problem
[D.Solution, D.Log] = Runfmincon( Obj, D0, [], 3 );
 
% Allocating Optimal Solution
x0.pdx1 = D.Solution.pdx1;
x0.pdx2 = D.Solution.pdx2;
x0.pdx3 = D.Solution.pdx3;

x0.ppx3 = D.Solution.ppx3;
x0.ppx4 = D.Solution.ppx4;

D.Surface = @(Pi, Fz, Gam, x0) ( x0.pdx1 + x0.pdx2.*dFz(Fz) ) .* ...
            ( 1 + x0.ppx3.*dPi(Pi) + x0.ppx4.*dPi(Pi).^2 ) .* ...
            ( 1 - x0.pdx3.*Gam.^2 ) .* Fz;

clear pdx1 pdx2 pdx3 ppx3 ppx4 D0

%% Fitting E_{x} Surface
% Initial Vector 
E0.pex1 = 0.1;
E0.pex2 = 0;
E0.pex3 = -0.1;

% Optimization Variables
pex1 = optimvar( 'pex1', 'Lowerbound',- 5   , 'Upperbound',  1   );
pex2 = optimvar( 'pex2', 'Lowerbound',- 2   , 'Upperbound',  1   );
pex3 = optimvar( 'pex3', 'Lowerbound',- 2   , 'Upperbound',- 0.01);

% Optimization Objective
Obj = fcn2optimexpr( @ErrorEx, pex1, pex2, pex3 );

% Optimization Constraint
Constr(1) = pex1 - pex2 + pex3 <= 1;
Constr(2) = pex1 - (pex2.^2)./(2*pex3) + (pex2.^2)./(4*pex3) <= 1;

% Solving Optimization Problem
[E.Solution, E.Log] = Runfmincon( Obj, E0, Constr, 3 );

% Allocating Solution
x0.pex1 = E.Solution.pex1;
x0.pex2 = E.Solution.pex2;
x0.pex3 = E.Solution.pex3;
x0.pex4 = 0;

E.Surface = @(Pi, Fz, Gam, x0) (x0.pex1 + x0.pex2.*dFz(Fz) + x0.pex3.*dFz(Fz).^2);

clear pex1 pex2 pex3 E0

%% Fitting K_{x \kappa} Surface
% Initial Vector 
K0.pkx1 = 50;
K0.pkx2 = 0;
K0.pkx3 = 0;

K0.ppx1 = 0;
K0.ppx2 = 0;

% Optimization Variables
pkx1 = optimvar( 'pkx1', 'Lowerbound', 0.1 );
pkx2 = optimvar( 'pkx2', 'Lowerbound', -5, 'Upperbound', 5 );
pkx3 = optimvar( 'pkx3', 'Lowerbound', -5, 'Upperbound', 5 );

ppx1 = optimvar( 'ppx1', 'Lowerbound', -5, 'Upperbound', 5 );
ppx2 = optimvar( 'ppx2', 'Lowerbound', -5, 'Upperbound', 5 );

% Optimization Objective
Obj = fcn2optimexpr( @ErrorKxk, pkx1, pkx2, pkx3, ppx1, ppx2 );

% Solving Optimization Problem
[K.Solution, K.Log] = Runfmincon( Obj, K0, [], 3 );

% Allocating Solution
x0.pkx1 = K.Solution.pkx1;
x0.pkx2 = K.Solution.pkx2;
x0.pkx3 = K.Solution.pkx3;

x0.ppx1 = K.Solution.ppx1;
x0.ppx2 = K.Solution.ppx2;

K.Surface = @(Pi, Fz, Gam, x0) Fz .* ( x0.pkx1 + x0.pkx2.*dFz(Fz) ) .* ...
    exp( x0.pkx3.*dFz(Fz) ) .* ( 1 + x0.ppx1.*dPi(Pi) + x0.ppx2.*dPi(Pi).^2 );

clear pkx1 pkx2 pkx3 ppx1 ppx2 K0

%% Fitting V_{x} Surface
% Initial Vector 
V0.pvx1 = 0;
V0.pvx2 = 0;

% Optimization Variables
pvx1 = optimvar( 'pvx1' );
pvx2 = optimvar( 'pvx2' );

% Optimization Objective
Obj = fcn2optimexpr( @ErrorVx, pvx1, pvx2 );

% Solving Optimization Problem
[V.Solution, V.Log] = Runfmincon( Obj, V0, [], 3 );

% Allocating Solution
x0.pvx1 = V.Solution.pvx1;
x0.pvx2 = V.Solution.pvx2;

V.Surface = @(Pi, Fz, Gam, x0) Fz .* ( ( x0.pvy1 + x0.pvy2.*dFz(Fz) ) );

clear pvx1 pvx2 V0

%% Fitting H_{x} Surface
% Initial Vector 
H0.phx1 = -0.2;
H0.phx2 = -0.2;

% Optimization Variables
phx1 = optimvar( 'phx1', 'Lowerbound', -5, 'Upperbound', 5 );
phx2 = optimvar( 'phx2', 'Lowerbound', -5, 'Upperbound', 5 );

% Optimization Objective
Obj = fcn2optimexpr( @ErrorHx, phx1, phx2 );

% Solving Optimization Problem
[H.Solution, H.Log] = Runfmincon( Obj, H0, [], 3 ); 

% Allocating Solution
x0.phx1 = H.Solution.phx1;
x0.phx2 = H.Solution.phx2;

H.Surface = @(Pi, Fz, Gam, x0) ( x0.phy1 + x0.phy2.*dFz(Fz) );

clear phx1 phx2 H0

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
            'MaxFunctionEvaluations', 10000, ...
            'MaxIterations', 10000, ...
            'Display', 'off', ...
            'PlotFcn', {@optimplotx, @optimplotfval}, ...
            'OutputFcn', @OptimLogging );
        
        % Solving Optimization Problem(s)
        for ii = 1 : n
            [Solution(ii), Feval(ii)] = solve( Prob, x0(ii), ...
                'solver', 'fmincon', 'options', Opts ); %#ok<AGROW>
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

    function RMSE = ErrorDx( pdx1, pdx2, pdx3, ppx3, ppx4 )
        DxSurface = (pdx1 + pdx2.*[Mesh.dFz]) .* ...
            (1 + ppx3.*[Mesh.dPi] + ppx4.*[Mesh.dPi].^2) .* ...
            (1 - pdx3.*[Mesh.Inclination].^2).*[Mesh.Load];
        
        RMSE = sqrt( mean( ([Nominal.D] - DxSurface).^2 ) );
    end

    function RMSE = ErrorEx( pex1, pex2, pex3 )
        ExSurface = (pex1 + pex2.*[Mesh.dFz] + pex3.*[Mesh.dFz].^2);
        
        RMSE = sqrt(mean( ([Nominal.E] - ExSurface).^2 ) );
    end

    function RMSE = ErrorKxk( pkx1, pkx2, pkx3, ppx1, ppx2 )
        KxkSurface = [Mesh.Load] .* ( pkx1 + pkx2.*[Mesh.dFz] ) .* ...
            exp( pkx3.*[Mesh.dFz] ) .* ( 1 + ppx1.*[Mesh.dPi] + ppx2.*[Mesh.dPi].^2 );
        
        RMSE = sqrt( mean( ([Nominal.K] - KxkSurface).^2 ) );
    end

    function RMSE = ErrorVx( pvx1, pvx2 )
        VxSurface = [Mesh.Load] .* ( pvx1 + pvx2.*[Mesh.dFz] );
        
        RMSE = sqrt( mean( ([Nominal.V] - VxSurface).^2 ) );
    end

    function RMSE = ErrorHx( phx1, phx2 )
        HxSurface = phx1 + phx2.*[Mesh.dFz];
        
        RMSE = sqrt( mean( ([Nominal.H] - HxSurface ).^2 ) );
    end
end