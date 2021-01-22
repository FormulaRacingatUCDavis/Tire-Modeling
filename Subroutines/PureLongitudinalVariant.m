function [ Variant, Tire ] = PureLongitudinalVariant( Raw, x0, Tire )

%% Optimization Variables
pcx1 = optimvar( 'pcx1', 'Lowerbound', 0 );

pdx1 = optimvar( 'pdx1', 'Lowerbound', 0 );
pdx2 = optimvar( 'pdx2', 'Upperbound', 0 );
pdx3 = optimvar( 'pdx3' );

pex1 = optimvar( 'pex1' );
pex2 = optimvar( 'pex2' );
pex3 = optimvar( 'pex3' );
pex4 = optimvar( 'pex4' );

pkx1 = optimvar( 'pkx1', 'Lowerbound', 0.1 );
pkx2 = optimvar( 'pkx2', 'Lowerbound', -5, 'Upperbound', 5 );
pkx3 = optimvar( 'pkx3', 'Lowerbound', -5, 'Upperbound', 5 );

phx1 = optimvar( 'phx1', 'Lowerbound', -5, 'Upperbound', 5 );
phx2 = optimvar( 'phx2', 'Lowerbound', -5, 'Upperbound', 5 );

pvx1 = optimvar( 'pvx1' );
pvx2 = optimvar( 'pvx2' );

ppx1 = optimvar( 'ppx1', 'Lowerbound', -5, 'Upperbound', 5 );
ppx2 = optimvar( 'ppx2', 'Lowerbound', -5, 'Upperbound', 5 );
ppx3 = optimvar( 'ppx3', 'Lowerbound', -5, 'Upperbound', 5 );
ppx4 = optimvar( 'ppx4', 'Lowerbound', -5, 'Upperbound', 5 );

%% Optimization Objective
Obj = fcn2optimexpr( @ErrorFyo, pcx1, ...
    pdx1, pdx2, pdx3, ...
    pex1, pex2, pex3, pex4, ...
    pkx1, pkx2, pkx3, ...
    phx1, phx2, ...
    pvx1, pvx2, ...
    ppx1, ppx2, ppx3, ppx4 );

%% Optimization Constraint
[Operating.dFz, Operating.Inclination] = meshgrid( (0:100:2500-Tire.Pacejka.Fzo)./Tire.Pacejka.Fzo, 0:5 );

Constr = optimineq( 2*numel(Operating.dFz) );

for i = 1 : numel(Operating.dFz)
    Constr(i) = ( pex1 + pex2.*Operating.dFz(i) + pex3.*Operating.dFz(i).^2 ) .* ...
        ( 1 - pex4.*Operating.Inclination(i) ) <= 0.95;
    
    Constr( i+numel(Operating.dFz) ) = ( pex1 + pex2.*Operating.dFz(i) + pex3.*Operating.dFz(i).^2 ) .* ...
        ( 1 + pex4.*Operating.Inclination(i) ) <= 0.95;
end

%% Solving Optimization Problem
[Variant.Solution, Variant.Log] = Runfmincon( Obj, x0, Constr, 3 );

%% Clearing Optimization Figure
delete( findobj( 'Type', 'figure', 'Name', 'Optimization PlotFcns' ) );

%% Allocating Solution
Tire.p.C.x = Variant.Solution.pcx1;

Tire.p.D.x(1) = Variant.Solution.pdx1;
Tire.p.D.x(2) = Variant.Solution.pdx2;
Tire.p.D.x(3) = Variant.Solution.pdx3;

Tire.p.E.x(1) = Variant.Solution.pex1;
Tire.p.E.x(2) = Variant.Solution.pex2;
Tire.p.E.x(3) = Variant.Solution.pex3;
Tire.p.E.x(4) = Variant.Solution.pex4;

Tire.p.K.x(1) = Variant.Solution.pkx1;
Tire.p.K.x(2) = Variant.Solution.pkx2;
Tire.p.K.x(3) = Variant.Solution.pkx3;

Tire.p.H.x(1) = Variant.Solution.phx1;
Tire.p.H.x(2) = Variant.Solution.phx2;

Tire.p.V.x(1) = Variant.Solution.pvx1;
Tire.p.V.x(2) = Variant.Solution.pvx2;

Tire.p.P.x(1) = Variant.Solution.ppx1;
Tire.p.P.x(2) = Variant.Solution.ppx2;
Tire.p.P.x(3) = Variant.Solution.ppx3;
Tire.p.P.x(4) = Variant.Solution.ppx4;

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

    function MeanSquareError = ErrorFyo( pcx1, ...
            pdx1, pdx2, pdx3, ...
            pex1, pex2, pex3, pex4, ...
            pkx1, pkx2, pkx3, ...
            phx1, phx2, ...
            pvx1, pvx2, ...
            ppx1, ppx2, ppx3, ppx4 )
        
        Cx = pcx1;
       
        Dx = (pdx1 + pdx2.*[Raw.dFz]) .* ...
            (1 + ppx3.*[Raw.dPi] + ppx4.*[Raw.dPi].^2) .* ...
            (1 - pdx3.*[Raw.Inclination].^2).*[Raw.Load];
        
        Ex = ( pex1 + pex2.*[Raw.dFz] + pex3.*[Raw.dFz].^2 ) .* ...
            ( 1 - pex4.*sign([Raw.Slip] ) );
        
        Kxk = [Raw.Load].*(pkx1 + pkx2.*[Raw.dFz]).*exp( pkx3.*[Raw.dFz] ).* ...
            (1 + ppx1.*[Raw.dPi] + ppx2.*[Raw.dPi].^2);
        
        Bx = Kxk ./ (Cx.*Dx);
        
        Vx = [Raw.Load].*(pvx1 + pvx2.*[Raw.dFz]);
        
        Hx = (phx1 + phx2.*[Raw.dFz]);
        
        Fxo = Dx.*sin( Cx.*atan( (1-Ex).*Bx.*([Raw.Slip] - Hx) + ...
            Ex.*atan(Bx.*([Raw.Slip] - Hx) ) ) ) + Vx;
        
        MeanSquareError = mean( ([Raw.Force] - Fxo).^2 );
    end
end