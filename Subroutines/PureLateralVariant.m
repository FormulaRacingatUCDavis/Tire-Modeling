function [ Variant, Tire ] = PureLateralVariant( Raw, x0, Tire )

%% Optimization Variables
pcy1 = optimvar( 'pcy1', 'Lowerbound', 0, 'Upperbound', 5 );

pdy1 = optimvar( 'pdy1', 'Lowerbound', 0, 'Upperbound', 10 );
pdy2 = optimvar( 'pdy2', 'Lowerbound', -5, 'Upperbound', 0 );
pdy3 = optimvar( 'pdy3', 'Lowerbound', -5, 'Upperbound', 5 );

pey1 = optimvar( 'pey1', 'Lowerbound', -5, 'Upperbound', 0.75 );
pey2 = optimvar( 'pey2', 'Lowerbound', -1.5, 'Upperbound', 1.5 );
pey3 = optimvar( 'pey3', 'Lowerbound', -1.5, 'Upperbound', 1.5 );
pey4 = optimvar( 'pey4', 'Lowerbound', -1.5, 'Upperbound', 1.5 );
pey5 = optimvar( 'pey5', 'Lowerbound', -1.5, 'Upperbound', 0 );

pky1 = optimvar( 'pky1', 'Lowerbound', 0.1, 'Upperbound', 25 );
pky2 = optimvar( 'pky2', 'Lowerbound', 0.1, 'Upperbound', 5 );
pky3 = optimvar( 'pky3', 'Lowerbound', -5, 'Upperbound', 5 );
pky4 = optimvar( 'pky4', 'Lowerbound', 0.05, 'Upperbound', 2.5 );
pky5 = optimvar( 'pky5', 'Lowerbound', -5, 'Upperbound', 5 );
pky6 = optimvar( 'pky6', 'Lowerbound', -5, 'Upperbound', 5 );
pky7 = optimvar( 'pky7', 'Lowerbound', -5, 'Upperbound', 5 );

phy1 = optimvar( 'phy1', 'Lowerbound', -5, 'Upperbound', 5 );
phy2 = optimvar( 'phy2', 'Lowerbound', -5, 'Upperbound', 5 );

pvy1 = optimvar( 'pvy1', 'Lowerbound', -5, 'Upperbound', 5 );
pvy2 = optimvar( 'pvy2', 'Lowerbound', -5, 'Upperbound', 5 );
pvy3 = optimvar( 'pvy3', 'Lowerbound', -5, 'Upperbound', 5 );
pvy4 = optimvar( 'pvy4', 'Lowerbound', -5, 'Upperbound', 5 );

ppy1 = optimvar( 'ppy1', 'Lowerbound', -5, 'Upperbound', 5 );
ppy2 = optimvar( 'ppy2', 'Lowerbound', -5, 'Upperbound', 5 );
ppy3 = optimvar( 'ppy3', 'Lowerbound', -5, 'Upperbound', 5 );
ppy4 = optimvar( 'ppy4', 'Lowerbound', -5, 'Upperbound', 5 );
ppy5 = optimvar( 'ppy5', 'Lowerbound', -5, 'Upperbound', 5 );

%% Optimization Objective
Obj = fcn2optimexpr( @ErrorFyo, pcy1, ...
    pdy1, pdy2, pdy3, ...
    pey1, pey2, pey3, pey4, pey5, ...
    pky1, pky2, pky3, pky4, pky5, pky6, pky7, ...
    phy1, phy2, ...
    pvy1, pvy2, pvy3, pvy4, ...
    ppy1, ppy2, ppy3, ppy4, ppy5 );

%% Optimization Constraint
Constr = optimineq( 4 );

Constr(1) = ( pey1 + pey2.*(0-Tire.Fzo)./Tire.Fzo ) .* ...
    ( 1 + pey5.*(pey4./(2.*pey5)).^2 - ...
    ( pey3 + pey4.*(pey4./(2.*pey5)) ) ) <= 0.99;

Constr(2) = ( pey1 + pey2.*(2500-Tire.Fzo)./Tire.Fzo ) .* ...
    ( 1 + pey5.*(pey4./(2.*pey5)).^2 - ...
    ( pey3 + pey4.*(pey4./(2.*pey5)) ) ) <= 0.99;

Constr(3) = ( pey1 + pey2.*(0-Tire.Fzo)./Tire.Fzo ) .* ...
    ( 1 + pey5.*(-pey4./(2.*pey5)).^2 + ...
    ( pey3 + pey4.*(-pey4./(2.*pey5)) ) ) <= 0.99;

Constr(4) = ( pey1 + pey2.*(2500-Tire.Fzo)./Tire.Fzo ) .* ...
    ( 1 + pey5.*(-pey4./(2.*pey5)).^2 + ...
    ( pey3 + pey4.*(-pey4./(2.*pey5)) ) ) <= 0.99;

%% Solving Optimization Problem
[Variant.Solution, Variant.Log] = Runfmincon( Obj, x0, Constr, 3 );

%% Clearing Optimization Figure
delete( findobj( 'Type', 'figure', 'Name', 'Optimization PlotFcns' ) );

%% Allocating Solution
Tire.p.C.y = Variant.Solution.pcy1;

Tire.p.D.y(1) = Variant.Solution.pdy1;
Tire.p.D.y(2) = Variant.Solution.pdy2;
Tire.p.D.y(3) = Variant.Solution.pdy3;

Tire.p.E.y(1) = Variant.Solution.pey1;
Tire.p.E.y(2) = Variant.Solution.pey2;
Tire.p.E.y(3) = Variant.Solution.pey3;
Tire.p.E.y(4) = Variant.Solution.pey4;
Tire.p.E.y(5) = Variant.Solution.pey5;

Tire.p.K.y(1) = Variant.Solution.pky1;
Tire.p.K.y(2) = Variant.Solution.pky2;
Tire.p.K.y(3) = Variant.Solution.pky3;
Tire.p.K.y(4) = Variant.Solution.pky4;
Tire.p.K.y(5) = Variant.Solution.pky5;
Tire.p.K.y(6) = Variant.Solution.pky6;
Tire.p.K.y(7) = Variant.Solution.pky7;

Tire.p.H.y(1) = Variant.Solution.phy1;
Tire.p.H.y(2) = Variant.Solution.phy2;

Tire.p.V.y(1) = Variant.Solution.pvy1;
Tire.p.V.y(2) = Variant.Solution.pvy2;
Tire.p.V.y(3) = Variant.Solution.pvy3;
Tire.p.V.y(4) = Variant.Solution.pvy4;

Tire.p.P.y(1) = Variant.Solution.ppy1;
Tire.p.P.y(2) = Variant.Solution.ppy2;
Tire.p.P.y(3) = Variant.Solution.ppy3;
Tire.p.P.y(4) = Variant.Solution.ppy4;
Tire.p.P.y(5) = Variant.Solution.ppy5;

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

    function MeanAbsoluteError = ErrorFyo( pcy1, ...
            pdy1, pdy2, pdy3, ...
            pey1, pey2, pey3, pey4, pey5, ...
            pky1, pky2, pky3, pky4, pky5, pky6, pky7, ...
            phy1, phy2, ...
            pvy1, pvy2, pvy3, pvy4, ...
            ppy1, ppy2, ppy3, ppy4, ppy5 )
        
        Cy = pcy1;
        
        Dy = (pdy1 + pdy2.*[Raw.dFz]) .* ...
            (1 + ppy3.*[Raw.dPi] + ppy4.*[Raw.dPi].^2) .* ...
            (1 - pdy3.*[Raw.Camber].^2).*[Raw.Load];
        
        Kya = pky1 .* Tire.Fzo .* ( 1 + ppy1.*[Raw.dPi] ) .* ...
            ( 1 - pky3.*abs([Raw.Camber]) ) .* sin( pky4 .* ...
            atan( ([Raw.Load]./Tire.Fzo) ./ ...
            ( ( pky2 + pky5.*[Raw.Camber].^2 ) .* ( 1 + ppy2.*[Raw.dPi] ) ) ) );
        
        Kyg0 = [Raw.Load].*(pky6 + pky7.*[Raw.dFz]) .* (1 + ppy5.*[Raw.dPi]);
        
        By = Kya ./ (Cy.*Dy);
        
        Vyg = [Raw.Load].*(pvy3 + pvy4.*[Raw.dFz]).*[Raw.Camber];
        
        Vy = [Raw.Load].*(pvy1 + pvy2.*[Raw.dFz]) + Vyg;
        
        Hy = (phy1 + phy2.*[Raw.dFz]) .* (Kyg0.*[Raw.Camber] - Vyg) ./ Kya;
        
        Ey = ( pey1 + pey2.*[Raw.dFz] ) .* ...
            ( 1 + pey5.*[Raw.Camber].^2 - ...
            ( pey3 + pey4.*[Raw.Camber] ).*sign([Raw.Slip] + Hy) );

        Fyo = Dy.*sin( Cy.*atan( (1-Ey).*By.*([Raw.Slip] + Hy) + ...
            Ey.*atan(By.*([Raw.Slip] + Hy) ) ) ) + Vy;
        
        MeanAbsoluteError = mean(abs([Raw.Force] - Fyo));
    end
end