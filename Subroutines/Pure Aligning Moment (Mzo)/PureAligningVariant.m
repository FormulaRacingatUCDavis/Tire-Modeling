function [ Variant, Tire ] = PureAligningVariant( Raw, x0, Tire )

%% Evaluate Fyo
Fyo = FyoEvaluation;

%% Optimization Variables
qbz1  = optimvar( 'qbz1' , 'Lowerbound',- 5   , 'Upperbound',  2    );
qbz2  = optimvar( 'qbz2' , 'Lowerbound',-Inf  , 'Upperbound',  2*x0.qbz2 );
qbz3  = optimvar( 'qbz3' , 'Lowerbound',-Inf  , 'Upperbound',  2*x0.qbz3 );
qbz5  = optimvar( 'qbz5' , 'Lowerbound',- 5   , 'Upperbound',  5    );
qbz6  = optimvar( 'qbz6' , 'Lowerbound',- 5   , 'Upperbound',  5    );
qbz10 = optimvar( 'qbz10', 'Lowerbound',- 5   , 'Upperbound',  5    );

qcz1  = optimvar( 'qcz1' , 'Lowerbound',  1   , 'Upperbound',  5    );

qdz1  = optimvar( 'qdz1' , 'Lowerbound',  0   , 'Upperbound', 10    );
qdz2  = optimvar( 'qdz2' , 'Lowerbound',- 1   , 'Upperbound',  0    );
qdz3  = optimvar( 'qdz3' , 'Lowerbound',- 5   , 'Upperbound',  5    );
qdz4  = optimvar( 'qdz4' , 'Lowerbound',- 5   , 'Upperbound',  5    );
qdz6  = optimvar( 'qdz6' , 'Lowerbound',- 5   , 'Upperbound',  5    );
qdz7  = optimvar( 'qdz7' , 'Lowerbound',- 1   , 'Upperbound',  0.1  );
qdz8  = optimvar( 'qdz8' , 'Lowerbound',- 5   , 'Upperbound',  5    );
qdz9  = optimvar( 'qdz9' , 'Lowerbound',- 1   , 'Upperbound',  0.1  );
qdz10 = optimvar( 'qdz10', 'Lowerbound',- 5   , 'Upperbound',  5    );
qdz11 = optimvar( 'qdz11', 'Lowerbound',- 1   , 'Upperbound',  0.1  );

qez1  = optimvar( 'qez1' , 'Lowerbound',- 5   , 'Upperbound',  0.8  );
qez2  = optimvar( 'qez2' , 'Lowerbound',- 5   , 'Upperbound',  0    );
qez3  = optimvar( 'qez3' , 'Lowerbound',- 0.5 , 'Upperbound',  0.1  );
qez4  = optimvar( 'qez4' , 'Lowerbound',- 1.5 , 'Upperbound',  1.5  );
qez5  = optimvar( 'qez5' , 'Lowerbound',- 1.5 , 'Upperbound',  0    );

qhz1  = optimvar( 'qhz1' , 'Lowerbound',- 5   , 'Upperbound',  5    );
qhz2  = optimvar( 'qhz2' , 'Lowerbound',- 5   , 'Upperbound',  5    );
qhz3  = optimvar( 'qhz3' , 'Lowerbound',- 5   , 'Upperbound',  5    );
qhz4  = optimvar( 'qhz4' , 'Lowerbound',- 5   , 'Upperbound',  5    );

ppz1  = optimvar( 'ppz1' , 'Lowerbound',- 5   , 'Upperbound',  5    );
ppz2  = optimvar( 'ppz2' , 'Lowerbound',- 5   , 'Upperbound',  5    );

%% Optimization Objective
Obj = fcn2optimexpr( @ErrorMzo, ...
    qbz1, qbz2, qbz3, qbz5, qbz6, qbz10, ...
    qcz1, ...
    qdz1, qdz2, qdz3, qdz4, qdz6, qdz7, qdz8, qdz9, qdz10, qdz11, ...
    qez1, qez2, qez3, qez4, qez5, ...
    qhz1, qhz2, qhz3, qhz4, ...
    ppz1, ppz2 );

%% Optimization Constraint
[dFz, Inclination] = meshgrid( ((0:50:2500)-Tire.Pacejka.Fzo)./Tire.Pacejka.Fzo, 0:0.1:5 );

Constr = optimineq( numel( dFz ) );

for i = 1 : numel( dFz )
    Constr(i) = fcn2optimexpr( @EtBound, ...
        qbz1, qbz2, qbz3, qbz5, qbz6, ...
        qcz1, ...
        qez1, qez2, qez3, qez4, qez5, ...
        dFz(i), Inclination(i) ) <= 0.99;
end

%% Solving Optimization Problem
[Variant.Solution, Variant.Log] = Runfmincon( Obj, x0, Constr, 3 );

%% Clearing Optimization Figure
delete( findobj( 'Type', 'figure', 'Name', 'Optimization PlotFcns' ) );

%% Allocating Solution
Tire.Pacejka.q.B.z(1)  = Variant.Solution.qbz1;
Tire.Pacejka.q.B.z(2)  = Variant.Solution.qbz2;
Tire.Pacejka.q.B.z(3)  = Variant.Solution.qbz3;
Tire.Pacejka.q.B.z(5)  = Variant.Solution.qbz5;
Tire.Pacejka.q.B.z(6)  = Variant.Solution.qbz6;
Tire.Pacejka.q.B.z(10) = Variant.Solution.qbz10;

Tire.Pacejka.q.C.z(1)  = Variant.Solution.qcz1;

Tire.Pacejka.q.D.z(1)  = Variant.Solution.qdz1;
Tire.Pacejka.q.D.z(2)  = Variant.Solution.qdz2;
Tire.Pacejka.q.D.z(3)  = Variant.Solution.qdz3;
Tire.Pacejka.q.D.z(4)  = Variant.Solution.qdz4;
Tire.Pacejka.q.D.z(6)  = Variant.Solution.qdz6;
Tire.Pacejka.q.D.z(7)  = Variant.Solution.qdz7;
Tire.Pacejka.q.D.z(8)  = Variant.Solution.qdz8;
Tire.Pacejka.q.D.z(9)  = Variant.Solution.qdz9;
Tire.Pacejka.q.D.z(10) = Variant.Solution.qdz10;
Tire.Pacejka.q.D.z(11) = Variant.Solution.qdz11;

Tire.Pacejka.q.E.z(1)  = Variant.Solution.qez1;
Tire.Pacejka.q.E.z(2)  = Variant.Solution.qez2;
Tire.Pacejka.q.E.z(3)  = Variant.Solution.qez3;
Tire.Pacejka.q.E.z(4)  = Variant.Solution.qez4;
Tire.Pacejka.q.E.z(5)  = Variant.Solution.qez5;

Tire.Pacejka.q.H.z(1)  = Variant.Solution.qhz1;
Tire.Pacejka.q.H.z(2)  = Variant.Solution.qhz2;
Tire.Pacejka.q.H.z(3)  = Variant.Solution.qhz3;
Tire.Pacejka.q.H.z(4)  = Variant.Solution.qhz4;

Tire.Pacejka.p.P.z(1)  = Variant.Solution.ppz1;
Tire.Pacejka.p.P.z(2)  = Variant.Solution.ppz2; 

%% Local Functions
    function Fyo = FyoEvaluation
        % Note Function is evaluated for null camber
        Fyo.Cy = Tire.Pacejka.p.C.y(1) .* ones( size( [Raw.Slip] ) );
        
        Fyo.Dy = (Tire.Pacejka.p.D.y(1) + Tire.Pacejka.p.D.y(2).*[Raw.dFz]) .* ...
            (1 + Tire.Pacejka.p.P.y(3).*[Raw.dPi] + Tire.Pacejka.p.P.y(4).*[Raw.dPi].^2) .* ...
            (1 - Tire.Pacejka.p.D.y(3).*0.^2).*[Raw.Load];
        
        Fyo.Kya = Tire.Pacejka.p.K.y(1) .* Tire.Pacejka.Fzo .* ( 1 + Tire.Pacejka.p.P.y(1).*[Raw.dPi] ) .* ...
            ( 1 - Tire.Pacejka.p.K.y(3).*abs(0) ) .* sin( Tire.Pacejka.p.K.y(4) .* ...
            atan( ([Raw.Load]./Tire.Pacejka.Fzo) ./ ...
            ( ( Tire.Pacejka.p.K.y(2) + Tire.Pacejka.p.K.y(5).*0.^2 ) .* ...
            ( 1 + Tire.Pacejka.p.P.y(2).*[Raw.dPi] ) ) ) );
        
        Fyo.Kyg0 = [Raw.Load].*(Tire.Pacejka.p.K.y(6) + Tire.Pacejka.p.K.y(7).*[Raw.dFz]) .* ...
            (1 + Tire.Pacejka.p.P.y(5).*[Raw.dPi]);
        
        Fyo.By = Fyo.Kya ./ (Fyo.Cy .* Fyo.Dy);
        
        Fyo.Vyg = [Raw.Load].*(Tire.Pacejka.p.V.y(3) + Tire.Pacejka.p.V.y(4).*[Raw.dFz]).*0;
        
        Fyo.Vy = [Raw.Load].*(Tire.Pacejka.p.V.y(1) + Tire.Pacejka.p.V.y(2).*[Raw.dFz]) + Fyo.Vyg;
        
        Fyo.Hy = (Tire.Pacejka.p.H.y(1) + Tire.Pacejka.p.H.y(2).*[Raw.dFz]) .* (Fyo.Kyg0.*0 - Fyo.Vyg) ./ Fyo.Kya;
        
        Fyo.Ey = ( Tire.Pacejka.p.E.y(1) + Tire.Pacejka.p.E.y(2).*[Raw.dFz] ) .* ...
            ( 1 + Tire.Pacejka.p.E.y(5).*0.^2 - ...
            ( Tire.Pacejka.p.E.y(3) + Tire.Pacejka.p.E.y(4).*0 ).*sign([Raw.Slip] + Fyo.Hy) );

        Fyo.Fyo = Fyo.Dy.*sin( Fyo.Cy.*atan( (1-Fyo.Ey).*Fyo.By.*([Raw.Slip] + Fyo.Hy) + ...
            Fyo.Ey.*atan(Fyo.By.*([Raw.Slip] + Fyo.Hy) ) ) ) + Fyo.Vy;
    end

    function [Solution, Log] = Runfmincon( Objective, x0, Constr, n )
        % Initializing Logs
        Log(n).x = [];
        Log(n).fval = [];
        
        % Creating Array of Initial Vectors (Latin Fyo.Hypercube Sampling)
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
            Prob = optimproblem( 'Objective', Objective );
        else
            Prob = optimproblem( 'Objective', Objective, 'Constraints', Constr );
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

    function Et = EtBound( qbz1, qbz2, qbz3, qbz5, qbz6, ...
            qcz1, ...
            qez1, qez2, qez3, qez4, qez5, ...
            dFz, Inclination )
        
        Et = (qez1 + qez2.*dFz + qez3.*dFz.^2) .* ...
            (1 + (qez4 + qez5.*Inclination).*(2/pi).* ...
            atan( (qbz1 + qbz2.*dFz + qbz3.*dFz^2) .* ...
            (1 + qbz5.*abs(Inclination) + qbz6.*Inclination^2) .* ...
            qcz1 .* deg2rad(45) ) );
    end

    function MeanAbsoluteError = ErrorMzo( ...
            qbz1, qbz2, qbz3, qbz5, qbz6, qbz10, ...
            qcz1, ...
            qdz1, qdz2, qdz3, qdz4, qdz6, qdz7, qdz8, qdz9, qdz10, qdz11, ...
            qez1, qez2, qez3, qez4, qez5, ...
            qhz1, qhz2, qhz3, qhz4, ...
            ppz1, ppz2 )
        
        Bt = (qbz1 + qbz2.*[Raw.dFz] + qbz3.*[Raw.dFz].^2) .* ...
            (1 + qbz5.*abs([Raw.Inclination]) + qbz6.*[Raw.Inclination].^2);
        
        Ct = qcz1;
        
        Dt = Tire.Pacejka.Ro .* ([Raw.Load]./Tire.Pacejka.Fzo) .* ...
            ( qdz1 + qdz2.*[Raw.dFz] ) .* ( 1 - ppz1.*[Raw.dPi] ) .* ...
            ( 1 + qdz3.*abs([Raw.Inclination]) + qdz4.*[Raw.Inclination].^2 );
        
        Ht = qhz1 + qhz2.*[Raw.dFz] + (qhz3 + qhz4.*[Raw.dFz]).*[Raw.Inclination];
        
        Et = (qez1 + qez2.*[Raw.dFz] + qez3.*[Raw.dFz].^2) .* ...
            (1 + (qez4 + qez5.*[Raw.Inclination]).*(2/pi).* ...
            atan( Bt .* Ct .* deg2rad([Raw.Slip]) ) );
        
        Br = qbz10 .* Fyo.By .* Fyo.Cy;
        
        Cr = 1;
        
        Dr = Tire.Pacejka.Ro .* [Raw.Load] .* ( ( qdz6 + qdz7.*[Raw.dFz] ) + ...
            ( ( qdz8 + qdz9.*[Raw.dFz] ) .* ( 1 + ppz2.*[Raw.dPi] ) + ...
            ( qdz10 + qdz11.*[Raw.dFz] ) .* abs([Raw.Inclination]) ).*[Raw.Inclination] );
        
        Hf = Fyo.Hy + Fyo.Vy ./ Fyo.Kya;
        
        % Evaluate Functions & Error
        t0 = Dt.*cos( Ct.*atan( (1-Et).*(Bt.*([Raw.Slip]+Ht) + ...
            Et.*atan( Bt.*([Raw.Slip]+Ht) ) ) ) ) .* cos( [Raw.Slip] );
        
        Mzro = Dr .* cos( Cr.*atan( Br.*([Raw.Slip]+Hf) ) ) .* cos( [Raw.Slip] );
        
        Mzo = -t0 .* Fyo.Fyo + Mzro;
        
        MeanAbsoluteError = mean(abs([Raw.Moment] - Mzo));
    end
end