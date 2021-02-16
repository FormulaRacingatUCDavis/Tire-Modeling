function [ Variant, Tire ] = PureAligningVariant( Tire, Raw, Nominal, Response )

%% Evaluate Fyo
[By, Cy, Kya, Hy, Vy, Fyo] = FyoEvaluation;

%% Optimization Variables
qbz1  = optimvar( 'qbz1' , 'Lowerbound',  0.1 , 'Upperbound', 10    );
qbz2  = optimvar( 'qbz2' , 'Lowerbound',- 1   , 'Upperbound',  1    );
qbz3  = optimvar( 'qbz3' , 'Lowerbound',- 1   , 'Upperbound',  1    );
qbz5  = optimvar( 'qbz5' , 'Lowerbound',- 0.5 , 'Upperbound',  0.5  );
qbz6  = optimvar( 'qbz6' , 'Lowerbound',- 0.1 , 'Upperbound',  0.1  );
qbz10 = optimvar( 'qbz10', 'Lowerbound',  0.01, 'Upperbound',  5    );

qcz1  = optimvar( 'qcz1' , 'Lowerbound',  0.9 , 'Upperbound',  1.1  );

qdz1  = optimvar( 'qdz1' , 'Lowerbound',  0.01, 'Upperbound',  5    );
qdz2  = optimvar( 'qdz2' , 'Lowerbound',- 0.75, 'Upperbound',  0    );
qdz3  = optimvar( 'qdz3' , 'Lowerbound',- 1   , 'Upperbound',  1    );
qdz4  = optimvar( 'qdz4' , 'Lowerbound',- 1   , 'Upperbound',  1    );

qdz6  = optimvar( 'qdz6' , 'Lowerbound',-Inf , 'Upperbound', Inf );
qdz7  = optimvar( 'qdz7' , 'Lowerbound',-Inf , 'Upperbound', Inf );
qdz8  = optimvar( 'qdz8' , 'Lowerbound',-Inf , 'Upperbound', Inf );
qdz9  = optimvar( 'qdz9' , 'Lowerbound',-Inf , 'Upperbound', Inf );
qdz10  = optimvar( 'qdz10' , 'Lowerbound',-Inf , 'Upperbound', Inf );
qdz11  = optimvar( 'qdz11' , 'Lowerbound',-Inf , 'Upperbound', Inf );
%qdz6  = optimvar( 'qdz6' , 'Lowerbound',-abs(x0.qdz6) , 'Upperbound', abs(x0.qdz6)  );
%qdz7  = optimvar( 'qdz7' , 'Lowerbound',-abs(x0.qdz7) , 'Upperbound', abs(x0.qdz7)  );
%qdz8  = optimvar( 'qdz8' , 'Lowerbound',-abs(x0.qdz8) , 'Upperbound', abs(x0.qdz8)  );
%qdz9  = optimvar( 'qdz9' , 'Lowerbound', abs(x0.qdz9) , 'Upperbound', abs(x0.qdz9)  );
%qdz10 = optimvar( 'qdz10', 'Lowerbound',-abs(x0.qdz10), 'Upperbound', abs(x0.qdz10) );
%qdz11 = optimvar( 'qdz11', 'Lowerbound',-abs(x0.qdz11), 'Upperbound', abs(x0.qdz11) );

qez1  = optimvar( 'qez1' , 'Lowerbound',-20   , 'Upperbound',- 1    );
qez2  = optimvar( 'qez2' , 'Lowerbound',-20   , 'Upperbound',  0    );
qez3  = optimvar( 'qez3' , 'Lowerbound',-25   , 'Upperbound',- 0.001);
qez4  = optimvar( 'qez4' , 'Lowerbound',- 1   , 'Upperbound',  1    );
qez5  = optimvar( 'qez5' , 'Lowerbound',- 1   , 'Upperbound',  1    );

qhz1  = optimvar( 'qhz1' , 'Lowerbound',- 0.03, 'Upperbound',  0.03 );
qhz2  = optimvar( 'qhz2' , 'Lowerbound',- 0.03, 'Upperbound',  0.03 );
qhz3  = optimvar( 'qhz3' , 'Lowerbound',- 0.01, 'Upperbound',  0.01 );
qhz4  = optimvar( 'qhz4' , 'Lowerbound',- 0.01, 'Upperbound',  0.01 );

ppz1  = optimvar( 'ppz1' , 'Lowerbound',- 1   , 'Upperbound',  1    );
ppz2  = optimvar( 'ppz2' , 'Lowerbound',- 1   , 'Upperbound',  1    );

%% Optimization Objective
Obj = fcn2optimexpr( @ErrorMzo, ...
    qbz1, qbz2, qbz3, qbz5, qbz6, qbz10, ...
    qcz1, ...
    qdz1, qdz2, qdz3, qdz4, qdz6, qdz7, qdz8, qdz9, qdz10, qdz11, ...
    qez1, qez2, qez3, qez4, qez5, ...
    qhz1, qhz2, qhz3, qhz4, ...
    ppz1, ppz2 );

%% Optimization Constraint
Constr(1) = fcn2optimexpr( @EtBound, qez1, qez2, qez3, qez4, qez5,  0, 0 ) <= 0;
Constr(2) = fcn2optimexpr( @EtBound, qez1, qez2, qez3, qez4, qez5, -1, 0 ) <= 0;
Constr(3) = fcn2optimexpr( @EtBound, qez1, qez2, qez3, qez4, qez5, -qez2./(2*qez3), 0 ) <= 0;

Constr(4) = fcn2optimexpr( @EtBound, qez1, qez2, qez3, qez4, qez5,  0, 5 ) <= 0;
Constr(5) = fcn2optimexpr( @EtBound, qez1, qez2, qez3, qez4, qez5, -1, 5 ) <= 0;
Constr(6) = fcn2optimexpr( @EtBound, qez1, qez2, qez3, qez4, qez5, -qez2./(2*qez3), 5 ) <= 0;

Constr(7) = fcn2optimexpr( @EtBound, qez1, qez2, qez3, qez4, qez5,  0, -5 ) <= 0;
Constr(8) = fcn2optimexpr( @EtBound, qez1, qez2, qez3, qez4, qez5, -1, -5 ) <= 0;
Constr(9) = fcn2optimexpr( @EtBound, qez1, qez2, qez3, qez4, qez5, -qez2./(2*qez3), -5 ) <= 0;

Fz = 0:100:Tire.Pacejka.Fzo;
Pi = [8 14];
Inc = -5:1:5;

[Fz,Pi,Inc] = meshgrid( Fz, Pi, Inc );

for i = 1:numel(Fz)
    Constr(i+9) = fcn2optimexpr( @DrBound, qdz6, qdz7, qdz8, qdz9, qdz10, qdz11, ppz2, ...
        Fz(i), Pi(i), Inc(i) ) <= max([Nominal.Dr]);
    Constr(i+9+numel(Fz)) = min([Nominal.Dr]) - ...
        fcn2optimexpr( @DrBound, qdz6, qdz7, qdz8, qdz9, qdz10, qdz11, ppz2, Fz(i), Pi(i), Inc(i) ) <= 0;
end

%% Setting Initialization
x0 = Response.x0;
x0.qdz6  = 0;
x0.qdz7  = 0;
x0.qdz8  = 0;
x0.qdz9  = 0;
x0.qdz10 = 0;
x0.qdz11 = 0;
x0.ppz2 = 0;

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
Tire.Pacejka.q.D.z(6)  = Variant.Solution.qdz6/10000;
Tire.Pacejka.q.D.z(7)  = Variant.Solution.qdz7/10000;
Tire.Pacejka.q.D.z(8)  = Variant.Solution.qdz8/10000;
Tire.Pacejka.q.D.z(9)  = Variant.Solution.qdz9/10000;
Tire.Pacejka.q.D.z(10) = Variant.Solution.qdz10/10000;
Tire.Pacejka.q.D.z(11) = Variant.Solution.qdz11/10000;

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
Tire.Pacejka.p.P.z(2)  = Variant.Solution.ppz2/10000; 

%% Local Functions
    function [By, Cy, Kya, Hy, Vy, Fyo] = FyoEvaluation
        % Fyo is Evaluated for Null Camber
        Cy = Tire.Pacejka.p.C.y(1) .* ones( size( [Raw.Slip] ) );
        
        Dy = (Tire.Pacejka.p.D.y(1) + Tire.Pacejka.p.D.y(2).*[Raw.dFz]) .* ...
            (1 + Tire.Pacejka.p.P.y(3).*[Raw.dPi] + Tire.Pacejka.p.P.y(4).*[Raw.dPi].^2) .* ...
            (1 - Tire.Pacejka.p.D.y(3).*0.^2).*[Raw.Load];
        
        Kya = Tire.Pacejka.p.K.y(1) .* Tire.Pacejka.Fzo .* ( 1 + Tire.Pacejka.p.P.y(1).*[Raw.dPi] ) .* ...
            ( 1 - Tire.Pacejka.p.K.y(3).*abs(0) ) .* sin( Tire.Pacejka.p.K.y(4) .* ...
            atan( ([Raw.Load]./Tire.Pacejka.Fzo) ./ ...
            ( ( Tire.Pacejka.p.K.y(2) + Tire.Pacejka.p.K.y(5).*0.^2 ) .* ...
            ( 1 + Tire.Pacejka.p.P.y(2).*[Raw.dPi] ) ) ) );
        
        Kyg0 = [Raw.Load].*(Tire.Pacejka.p.K.y(6) + Tire.Pacejka.p.K.y(7).*[Raw.dFz]) .* ...
            (1 + Tire.Pacejka.p.P.y(5).*[Raw.dPi]);
        
        By = Kya ./ (Cy .* Dy);
        
        Vyg = [Raw.Load].*(Tire.Pacejka.p.V.y(3) + Tire.Pacejka.p.V.y(4).*[Raw.dFz]).*0;
        
        Vy = [Raw.Load].*(Tire.Pacejka.p.V.y(1) + Tire.Pacejka.p.V.y(2).*[Raw.dFz]) + Vyg;
        
        Hy = (Tire.Pacejka.p.H.y(1) + Tire.Pacejka.p.H.y(2).*[Raw.dFz]) .* (Kyg0.*0 - Vyg) ./ Kya;
        
        Ey = ( Tire.Pacejka.p.E.y(1) + Tire.Pacejka.p.E.y(2).*[Raw.dFz] ) .* ...
            ( 1 + Tire.Pacejka.p.E.y(5).*0.^2 - ...
            ( Tire.Pacejka.p.E.y(3) + Tire.Pacejka.p.E.y(4).*0 ).*sign([Raw.Slip] + Hy) );

        Fyo = Dy.*sin( Cy.*atan( (1-Ey).*By.*([Raw.Slip] + Hy) + ...
            Ey.*atan(By.*([Raw.Slip] + Hy) ) ) ) + Vy;
        
        % Other Required Parameters Evaluated with Camber
        Dy = (Tire.Pacejka.p.D.y(1) + Tire.Pacejka.p.D.y(2).*[Raw.dFz]) .* ...
            (1 + Tire.Pacejka.p.P.y(3).*[Raw.dPi] + Tire.Pacejka.p.P.y(4).*[Raw.dPi].^2) .* ...
            (1 - Tire.Pacejka.p.D.y(3).*[Raw.Inclination].^2).*[Raw.Load];
        
        Kya = Tire.Pacejka.p.K.y(1) .* Tire.Pacejka.Fzo .* ( 1 + Tire.Pacejka.p.P.y(1).*[Raw.dPi] ) .* ...
            ( 1 - Tire.Pacejka.p.K.y(3).*abs([Raw.Inclination]) ) .* sin( Tire.Pacejka.p.K.y(4) .* ...
            atan( ([Raw.Load]./Tire.Pacejka.Fzo) ./ ...
            ( ( Tire.Pacejka.p.K.y(2) + Tire.Pacejka.p.K.y(5).*[Raw.Inclination].^2 ) .* ...
            ( 1 + Tire.Pacejka.p.P.y(2).*[Raw.dPi] ) ) ) );
        
        Kyg0 = [Raw.Load].*(Tire.Pacejka.p.K.y(6) + Tire.Pacejka.p.K.y(7).*[Raw.dFz]) .* ...
            (1 + Tire.Pacejka.p.P.y(5).*[Raw.dPi]);
        
        By = Kya ./ (Cy .* Dy);
        
        Vyg = [Raw.Load].*(Tire.Pacejka.p.V.y(3) + Tire.Pacejka.p.V.y(4).*[Raw.dFz]).*[Raw.Inclination];
        
        Vy = [Raw.Load].*(Tire.Pacejka.p.V.y(1) + Tire.Pacejka.p.V.y(2).*[Raw.dFz]) + Vyg;
        
        Hy = (Tire.Pacejka.p.H.y(1) + Tire.Pacejka.p.H.y(2).*[Raw.dFz]) .* ...
            (Kyg0.*[Raw.Inclination] - Vyg) ./ Kya;
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
            'MaxFunctionEvaluations', 10000000, ...
            'MaxIterations', 10000, ...
            'Display', 'off', ...
            'PlotFcn', {@optimplotx, @optimplotfval, @optimplotconstrviolation}, ...
            'OutputFcn', @OptimLogging );
        
        % Solving Optimization Problem(s)
        for ii = 1 : n
            % try
                [Solution(ii), Feval(ii)] = solve( Prob, x0(ii), ...
                    'solver', 'fmincon', 'options', Opts ); %#ok<AGROW>
            % catch
            %     continue
            % end
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

    function Et = EtBound( qez1, qez2, qez3, qez4, qez5, dFz, Inc )
        if (dFz > 0) | (dFz < -1)
            Et = 0;
        else
            Et = (qez1 + qez2.*dFz + qez3.*dFz.^2) .* ...
                (1 + (qez4 + qez5.*Inc ) ); 
        end
    end

    function Dr = DrBound( qdz6, qdz7, qdz8, qdz9, qdz10, qdz11, ppz2, Fz, Pi, Inc)
        dFz = (Fz - Tire.Pacejka.Fzo) ./ Tire.Pacejka.Fzo;
        dPi = (Pi - Tire.Pacejka.Pio) ./ Tire.Pacejka.Pio;
        
        Dr = Tire.Pacejka.Ro .* Fz .* ( ( qdz6 + qdz7.*dFz ) + ...
            ( ( qdz8 + qdz9.*dFz ) .* ( 1 + ppz2.*dPi ) + ...
            ( qdz10 + qdz11.*dFz ) .* abs(Inc) ).*Inc );
    end

    function RMSE = ErrorMzo( ...
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
        
        Br = qbz10 .* By .* Cy;
        
        Cr = 1;
        
        Dr = Tire.Pacejka.Ro .* [Raw.Load] .* ( ( qdz6/10000 + qdz7/10000.*[Raw.dFz] ) + ...
            ( ( qdz8/10000 + qdz9/10000.*[Raw.dFz] ) .* ( 1 + ppz2/10000.*[Raw.dPi] ) + ...
            ( qdz10/10000 + qdz11/10000.*[Raw.dFz] ) .* abs([Raw.Inclination]) ).*[Raw.Inclination] ) .* cos([Raw.Slip]);
        
        Hf = Hy + Vy ./ Kya;
        
        % Evaluate Functions & Error
        t0 = Dt.*cos( Ct.*atan( (1-Et).*(Bt.*([Raw.Slip]+Ht)) + ...
            Et.*atan( Bt.*([Raw.Slip]+Ht) ) ) ) .* cos([Raw.Slip]);
        
        Mzro = Dr .* cos( Cr.*atan( Br.*([Raw.Slip]+Hf) ) ) .* cos([Raw.Slip]);
        
        Mzo = -t0 .* Fyo + Mzro;
        
        RMSE = sqrt( mean( ([Raw.Moment] - Mzo).^2 ) );
    end
end