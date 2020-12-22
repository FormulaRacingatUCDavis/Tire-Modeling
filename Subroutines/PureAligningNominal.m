function [Nominal] = PureAligningNominal( Mesh, Raw, Tire )

%% Evaluating Pure Lateral Force Fit (Fyo)
Fyo = FyoEvaluation( Mesh, Raw, Tire );

%% Defining Initial Point
x0.Bt = 2.5;
x0.Ct = 0.03;
x0.Dt = -5;
x0.Et = 0;
x0.Ht = -0.15;

x0.qbz10 = 0.1;
x0.Dr = 0.1;

%% Defining Optimization Variables
Bt = optimvar( 'Bt', 'LowerBound', 0.001, 'UpperBound', 50);
Ct = optimvar( 'Ct', 'LowerBound', 0.001, 'UpperBound', 0.8);
Dt = optimvar( 'Dt', 'LowerBound', -2* max(abs(AligningMoment)), 'UpperBound',-0.1 );
Et = optimvar( 'Et', 'Lowerbound', -3, 'UpperBound', 0.9);
Ht = optimvar( 'Ht', 'Lowerbound', -2.5, 'UpperBound', 0);

qbz10 = optimvar( 'qbz10', 'LowerBound', 0 );
Cr = 1;
Dr = optimvar( 'Dr', 'LowerBound', 0, 'UpperBound', max(abs(Raw.Moment))/10 );
Hf = Fyo.Hy + Fyo.Vy / Fyo.Kya;

%% Defining Optimization Objective
Obj = fcn2optimexpr( @ErrorMzoNom, Bt, Ct, Dt, Et, Ht, qbz10, Cr, Dr, Hf );

%% Defining Optimization Problem
Prob = optimproblem( 'Objective', Obj );

%% Defining Optimization Options
Opts = optimoptions( 'fmincon', ...
    'MaxFunctionEvaluations', 10000, ...
    'MaxIterations', 10000, ...
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

Nominal.Residual = FyoNomResidual( Nominal.C, Nominal.D, Nominal.E, ...
    Nominal.K, Nominal.H, Nominal.V, [Raw.Slip; Raw.Force] );

%% Local Functions
    function Fyo = FyoEvaluation( Mesh, Raw, Tire ) 
        Fyo.Cy = Tire.p.C.y(1);
        
        Fyo.Dy = (Tire.p.D.y(1) + Tire.p.D.y(2).*Mesh.dFz) .* ...
            (1 + Tire.p.P.y(3).*Mesh.dPi + Tire.p.P.y(4).*Mesh.dPi.^2) .* ...
            (1 - Tire.p.D.y(3).*0.^2).*Mesh.Load;
        
        Fyo.Ey = ( Tire.p.E.y(1) + Tire.p.E.y(2).*Mesh.dFz ) .* ...
            ( 1 + Tire.p.E.y(5).*0.^2 - ...
            ( Tire.p.E.y(3) + Tire.p.E.y(4).*0 ).*sign(Raw.Slip) );
        
        Fyo.Kya = Tire.p.K.y(1) .* Tire.Fzo .* ( 1 + Tire.p.P.y(1).*Mesh.dPi ) .* ...
            ( 1 - Tire.p.K.y(3).*abs(0) ) .* sin( Tire.p.K.y(4) .* ...
            atan( (Mesh.Load./Tire.Fzo) ./ ...
            ( ( Tire.p.K.y(2) + Tire.p.K.y(5).*0.^2 ) .* ( 1 + Tire.p.P.y(2).*Mesh.dPi ) ) ) );
        
        Fyo.Kyg0 = Mesh.Load.*(Tire.p.K.y(6) + Tire.p.K.y(7).*Mesh.dFz ) .* (1 + Tire.p.P.y(5).*Mesh.dPi );
        
        Fyo.By = Fyo.Kya ./ ( Fyo.Cy.*Fyo.Dy );
        
        Fyo.Vyg = Mesh.Load.*(Tire.p.V.y(3) + Tire.p.V.y(4).*Mesh.dFz ).*0;
        
        Fyo.Vy = Mesh.Load.*(Tire.p.V.y(1) + Tire.p.V.y(2).*Mesh.dFz ) + Fyo.Vyg;
        
        Fyo.Hy = (Tire.p.H.y(1) + Tire.p.H.y(2).*Mesh.dFz ) .* ...
            (Fyo.Kyg0.*0 - Fyo.Vyg ) ./ Fyo.Kya;
        
        Fyo.Feval = Fyo.Dy.* sin( Fyo.Cy .* atan( (1-Fyo.Ey) .* ...
            Fyo.By.*(Raw.Slip - Fyo.Hy ) + Fyo.Ey.*atan( ...
            Fyo.By.*(Raw.Slip - Fyo.Hy ) ) ) ) + Fyo.Vy;
    end

    function MeanSquareError = ErrorMzoNom( Bt, Ct, Dt, Et, Ht, qbz10, Cr, Dr, Hf )
        % Magic Formula Evaluation
        MzoNom = -(Dt .* cos( Ct .* atan( (1-Et) .* ...
            (Bt.*(Raw.Slip-Ht)) + Et.*atan(Bt.*(Raw.Slip-Ht)) ) ) ) .* ...
            (Fyo.Dy .* sin( Fyo.Cy .* atan( (1-Fyo.Ey) .* ...
            (Fyo.By.*(Raw.Slip-Fyo.Hy)) + Fyo.Ey .* atan(Fyo.By.*(Raw.Slip-Fyo.Hy) ) ) ) + Fyo.Vy) + ...
            (Dr .* cos( Cr .* atan( qbz10 .* Fyo.By .* Fyo.Cy .* ...
            (Raw.Slip-Hf) ) ).*cos(Raw.Slip) );
        
        % Mean Square Error Evaluation
        MeanSquareError = mean( ( Raw.Moment - MzoNom ).^2 );
    end

    function Residual = FyoNomResidual( C, D, E, K, H, V, Data )
        % Magic Formula Evaluation
        B = K ./ ( C.*D + 0.1);
        
        Fit = ( D.*sin( C.*atan( ...
            ( 1 - E ).*( B.*( Data(1,:) - H ) ) + ...
            E.*atan( B .* ( Data(1,:) - H ) ) ) ) ) + V;
        
        % Mean Square Error Evaluation
        Residual = Data(2,:) - Fit;
    end
end


