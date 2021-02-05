function [Nominal] = PureLongitudinalNominal( Raw )

%% Defining Initial Point
[x0, lb, ub] = FxoInitialBounds( Raw.Slip, Raw.Force );

%% Defining Optimization Variables
C = optimvar( 'C', 'LowerBound', lb.C, 'UpperBound', ub.C );
D = optimvar( 'D', 'LowerBound', lb.D, 'UpperBound', ub.D );
E = optimvar( 'E', 'LowerBound', lb.E, 'UpperBound', ub.E );
K = optimvar( 'K', 'LowerBound', lb.K, 'UpperBound', ub.K );
H = optimvar( 'H', 'LowerBound', lb.H, 'UpperBound', ub.H );
V = optimvar( 'V', 'LowerBound', lb.V, 'UpperBound', ub.V );

%% Defining Optimization Objective
Obj = fcn2optimexpr( @ErrorFxoNom, C, D, E, K, H, V, [Raw.Slip; Raw.Force] );

%% Defining Optimization Problem
Prob = optimproblem( 'Objective', Obj );

%% Defining Optimization Options
Opts = optimoptions( 'fmincon', ...
    'MaxFunctionEvaluations', 2000, ...
    'MaxIterations', 2000, ...
    'Display', 'off' );
    
%% Solving Optimization Problem
Solution = solve( Prob, x0, 'solver', 'fmincon', 'options', Opts );

%% Allocating Solution
Nominal.C0 = x0.C;
Nominal.D0 = x0.D;
Nominal.E0 = x0.E;
Nominal.K0 = x0.K;
Nominal.H0 = x0.H;
Nominal.V0 = x0.V;

Nominal.C = Solution.C;
Nominal.D = Solution.D;
Nominal.E = Solution.E;
Nominal.K = Solution.K;
Nominal.H = Solution.H;
Nominal.V = Solution.V;

Nominal.Residual = FxoNomResidual( Nominal.C, Nominal.D, Nominal.E, ...
    Nominal.K, Nominal.H, Nominal.V, [Raw.Slip; Raw.Force] );

%% Local Functions
    function [x0, lb, ub] = FxoInitialBounds( SlipRatio, LongitudinalForce )
        [ ~, MaxIdx ] = maxk( LongitudinalForce, 5 );
        [~, MinIdx ] = mink( LongitudinalForce, 5 );
        
        MaxIdx = round( mean( MaxIdx ) );
        MinIdx = round( mean( MinIdx ) );

        LinearFit = fitlm( SlipRatio( (SlipRatio > -0.03) & (SlipRatio < 0.03) ), ...
            LongitudinalForce( (SlipRatio > -0.03) & (SlipRatio < 0.03) ) );
                
        x0.V = mean( [max(LongitudinalForce) min(LongitudinalForce)] );
        x0.K = LinearFit.Coefficients.Estimate(2);
        x0.H = ( x0.V - LinearFit.Coefficients.Estimate(1) ) ./ x0.K;
        
        xm0 = mean( abs( [ SlipRatio(MaxIdx) - x0.H, SlipRatio(MinIdx) - x0.H ] ) );
        
        x0.D = abs( max(LongitudinalForce) ) - x0.V;
        x0.C = 1 + ( 1 - 2*acos(0.8)/pi );
        
        B0 = x0.K / ( x0.C*x0.D );

        x0.E = ( B0*xm0 - tand( pi/(2*x0.C) ) ) / ( B0*xm0 - atand(B0*xm0) );
        
        %% Bounding
        lb.C = 0.5;
        lb.D = 0;
        lb.E = -Inf;
        lb.K = 0;
        lb.H = -Inf;
        lb.V = -Inf;
        
        ub.C = 1.5;
        ub.D = 1.5*max( abs( LongitudinalForce ) );
        ub.E = 0.95;
        ub.K = Inf;
        ub.H = Inf;
        ub.V = Inf;
    end
    
    function MeanSquareError = ErrorFxoNom( C, D, E, K, H, V, Data )
        % Magic Formula Evaluation
        B = K ./ ( C.*D + 0.1);
        
        Fit = ( D.*sin( C.*atan( ...
            ( 1 - E ).*( B.*( Data(1,:) + H ) ) + ...
            E.*atan( B .* ( Data(1,:) + H ) ) ) ) ) + V;
        
        % Mean Square Error Evaluation
        MeanSquareError = mean( ( Data(2,:) - Fit ).^2 );
    end

    function Residual = FxoNomResidual( C, D, E, K, H, V, Data )
        % Magic Formula Evaluation
        B = K ./ ( C.*D + 0.1);
        
        Fit = ( D.*sin( C.*atan( ...
            ( 1 - E ).*( B.*( Data(1,:) + H ) ) + ...
            E.*atan( B .* ( Data(1,:) + H ) ) ) ) ) + V;
        
        % Mean Square Error Evaluation
        Residual = Data(2,:) - Fit;
    end
end