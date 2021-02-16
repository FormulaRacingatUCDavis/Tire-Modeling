 function [Nominal] = PureAligningNominal( Tire, Raw, Mesh )
%% Evaluating Pure Lateral Force Fit (Fyo)
[Fyo, By, Cy, Kya, Hy, Vy] = FyoEval( Tire );

%% Optimization Variables
Bt = optimvar( 'Bt', 'Lowerbound',  0   , 'Upperbound', 10    );
Ct = optimvar( 'Ct', 'Lowerbound',  0.9 , 'Upperbound',  1.1  );
Dt = optimvar( 'Dt', 'Lowerbound', 0.005, 'Upperbound', Inf   );
Et = optimvar( 'Et', 'Lowerbound',- 25  , 'Upperbound',  0.95 );
Ht = optimvar( 'Ht', 'Lowerbound',- 0.04, 'Upperbound',  0.04 );

qbz10 = optimvar( 'qbz10', 'Lowerbound', 0.01, 'Upperbound', 2 );
Dr = optimvar( 'Dr', 'Lowerbound',-max(abs(Raw.Moment))/20, 'Upperbound', max(abs(Raw.Moment))/20 );

%% Optimization Objective
Obj = fcn2optimexpr( @ErrorMzo, Bt, Ct, Dt, Et, Ht, qbz10, Dr );

%% Solving Optimization Problem
x0.Bt = 1;
x0.Ct = 1;
x0.Dt = 0.05;
x0.Et = -2;
x0.Ht = 0;

x0.qbz10 = 1;
x0.Dr = 0;

[Fit.Solution, Fit.Log] = Runfmincon( Obj, x0, [], 1 );

%% Clearing Optimization Figure
delete( findobj( 'Type', 'figure', 'Name', 'Optimization PlotFcns' ) );

%% Allocating Solution
Nominal.Bt = Fit.Solution.Bt;
Nominal.Ct = Fit.Solution.Ct;
Nominal.Dt = Fit.Solution.Dt;
Nominal.Et = Fit.Solution.Et;
Nominal.Ht = Fit.Solution.Ht;

Nominal.qbz10 = Fit.Solution.qbz10;
Nominal.Cr = 1;
Nominal.Dr = Fit.Solution.Dr;
Nominal.Hf = Hy( Mesh.Pressure, Mesh.Load, Mesh.Inclination) + ...
             Vy( Mesh.Load, Mesh.Inclination) ./ Kya( Mesh.Pressure, Mesh.Load, Mesh.Inclination);
         
Nominal.Residual = ( Raw.Moment - MzoEval( Raw.Slip ) )';

%{
figure( 'Name', 'Aligning Moment Test', 'NumberTitle', 'off' )
yyaxis left
scatter( rad2deg(Raw.Slip), Raw.Moment, 'k.'); hold on;
fplot( @(Slip) MzoEval( deg2rad(Slip) ), [-15 15], 'g-.')
scatter( rad2deg(Raw.Slip), Nominal.Residual, 'r.' )
yyaxis right
fplot( @(Slip) Fyo(Mesh.Pressure, Mesh.Load, Mesh.Inclination, deg2rad(Slip)), [-15 15], 'b' )
pause(0.25)
delete( gcf )
%}

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

    function Residual = ErrorMzo( Bt, Ct, Dt, Et, Ht, qbz10, Dr )
        t0 = Dt.*cos( Ct.*atan( (1 - Et) .* Bt.*(Raw.Slip + Ht) + ...
            Et.*atan( Bt.*(Raw.Slip + Ht) ) ) ) .* cos(Raw.Slip);

        Hf = Hy( Raw.Pressure, Raw.Load, Raw.Inclination) + ...
            Vy( Raw.Load, Raw.Inclination) ./ Kya( Raw.Pressure, Raw.Load, Raw.Inclination);

        Mzr0 = Dr.*cos( atan( qbz10.*By( Raw.Pressure, Raw.Load, Raw.Inclination).*Cy .* ...
            (Raw.Slip + Hf) ) ) .* cos(Raw.Slip) .* cos(Raw.Slip);

        Mz0 = -t0 .* Fyo( Raw.Pressure, Raw.Load, 0, Raw.Slip ) + Mzr0;

        Residual = sqrt( mean( (Raw.Moment - Mz0).^2 ) );
    end
    
    function Mz0 = MzoEval( Slip )
        t0 = Nominal.Dt.*cos( Nominal.Ct.*atan( (1 - Nominal.Et) .* Nominal.Bt.*(Slip + Nominal.Ht) + ...
            Nominal.Et.*atan( Nominal.Bt.*(Slip + Nominal.Ht) ) ) ) .* cos(Slip);

        Mzr0 = Nominal.Dr.*cos( atan( Nominal.qbz10.*By( Mesh.Pressure, Mesh.Load, Mesh.Inclination).*Cy .* ...
            (Slip + Nominal.Hf) ) ) .* cos(Slip);

        Mz0 = -t0 .* Fyo( Mesh.Pressure, Mesh.Load, 0, Slip ) + Mzr0;
    end

    function [Fyo, By, Cy, Kya, Hy, Vy] = FyoEval( Tire )
        % Operating Condition Functions
        dFz = @(Fz) (Fz - Tire.Pacejka.Fzo) ./ Tire.Pacejka.Fzo;
        dPi = @(Pi) (Pi - Tire.Pacejka.Pio) ./ Tire.Pacejka.Pio;
        
        % Lateral Force Evaluation
        Cy = Tire.Pacejka.p.C.y(1);
        
        Dy = @(Pi, Fz, Gam) (Tire.Pacejka.p.D.y(1) + Tire.Pacejka.p.D.y(2).*dFz(Fz)) .* ...
            (1 + Tire.Pacejka.p.P.y(3).*dPi(Pi) + Tire.Pacejka.p.P.y(4).*dPi(Pi).^2) .* ...
            (1 - Tire.Pacejka.p.D.y(3).*Gam.^2).*Fz;
        
        Kya = @(Pi, Fz, Gam) Tire.Pacejka.p.K.y(1) .* Tire.Pacejka.Fzo .* ( 1 + Tire.Pacejka.p.P.y(1).*dPi(Pi) ) .* ...
            ( 1 - Tire.Pacejka.p.K.y(3).*abs(Gam) ) .* sin( Tire.Pacejka.p.K.y(4) .* ...
            atan( (Fz./Tire.Pacejka.Fzo) ./ ...
            ( ( Tire.Pacejka.p.K.y(2) + Tire.Pacejka.p.K.y(5).*Gam.^2 ) .* ( 1 + Tire.Pacejka.p.P.y(2).*dPi(Pi) ) ) ) );
        
        Kyg0 = @(Pi, Fz) Fz.*(Tire.Pacejka.p.K.y(6) + Tire.Pacejka.p.K.y(7).*dFz(Fz) ) .* (1 + Tire.Pacejka.p.P.y(5).*dPi(Pi) );
        
        By = @(Pi, Fz, Gam) Kya(Pi, Fz, 0) ./ ( Cy.*Dy(Pi, Fz, 0) );
        
        Vyg = @(Fz, Gam) Fz.*(Tire.Pacejka.p.V.y(3) + Tire.Pacejka.p.V.y(4).*dFz(Fz) ).*Gam;
        
        Vy = @(Fz, Gam) Fz.*(Tire.Pacejka.p.V.y(1) + Tire.Pacejka.p.V.y(2).*dFz(Fz) ) + Vyg(Fz, Gam);
        
        Hy = @(Pi, Fz, Gam) (Tire.Pacejka.p.H.y(1) + Tire.Pacejka.p.H.y(2).*dFz(Fz) ) .* ...
            (Kyg0(Pi, Fz).*Gam - Vyg(Fz, Gam) ) ./ Kya(Pi, Fz, Gam);
        
        Ey = @(Fz, Gam, Slip, Hy) ( Tire.Pacejka.p.E.y(1) + Tire.Pacejka.p.E.y(2).*dFz(Fz) ) .* ...
            ( 1 + Tire.Pacejka.p.E.y(5).*Gam.^2 - ...
            ( Tire.Pacejka.p.E.y(3) + Tire.Pacejka.p.E.y(4).*Gam ).*sign(Slip + Hy) );

        Fyo = @(Pi, Fz, Gam, Slip) Dy(Pi, Fz, Gam) .* ...
            sin( Cy .* atan( (1-Ey(Fz, Gam, Slip, Hy(Pi, Fz, Gam) )) .* ...
            By(Pi, Fz, Gam).*(Slip + Hy(Pi, Fz, Gam) ) + ...
            Ey(Fz, Gam, Slip, Hy(Pi, Fz, Gam) ).*atan( ...
            By(Pi, Fz, Gam).*(Slip + Hy(Pi, Fz, Gam) ) ) ) ) + Vy(Fz, Gam);
    end
end
