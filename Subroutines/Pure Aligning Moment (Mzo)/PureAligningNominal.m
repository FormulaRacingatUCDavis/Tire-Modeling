 function [Nominal] = PureAligningNominal( Mesh, Raw, Tire )
%% Evaluating Pure Lateral Force Fit (Fyo)
Fyo = FyoEvaluation( Mesh, Tire );

%% Preprocess Raw Data
[Data.Slip, Data.Moment] = prepareCurveData( Raw.Slip, Raw.Moment );

% Sort for Monotonicity
[Data.Slip, SortIdx] = sort( Data.Slip );
Data.Moment = Data.Moment( SortIdx );

%% Fit Type
Fit.Type = fittype( ['-(Dt*cos( Ct*atan((1-Et)*(Bt*(Slip+Ht)) + Et*atan(Bt*(Slip+Ht))))*cos(Slip)) * ', ... % Pneumatic Trail Function
    ... % Negative \alpha_{y} Fyo
    '( (Slip <= ', num2str(Fyo.Hy), ')*(', num2str(Fyo.Dy), '*sin(', num2str(Fyo.Cy), '*atan((1-', num2str(Fyo.Ey(1)), ...
    ')*(', num2str(Fyo.By), '*(Slip+', num2str(Fyo.Hy), '))+', num2str(Fyo.Ey(1)), ...
    '*atan(', num2str(Fyo.By), '*(Slip+', num2str(Fyo.Hy), ')))) +', num2str(Fyo.Vy), ') + ' , ...
    ... % Positive \alpha_{y} Fyo
    '(Slip > ', num2str(Fyo.Hy), ')*(', num2str(Fyo.Dy), '*sin(', num2str(Fyo.Cy), '*atan((1-', num2str(Fyo.Ey(2)), ...
    ')*(', num2str(Fyo.By), '*(Slip+', num2str(Fyo.Hy), '))+', num2str(Fyo.Ey(2)), ...
    '*atan(', num2str(Fyo.By), '*(Slip+', num2str(Fyo.Hy), ')))) +', num2str(Fyo.Vy), ') ) + ', ...
    ... % Residual Torque Function
    '(Dr*cos(atan(qbz10*', num2str(Fyo.By*Fyo.Cy), '*(Slip+', num2str(Fyo.Hy+Fyo.Vy/Fyo.Kya), ')))*cos(Slip))'], ...
    ... % Fit Options
    'independent', 'Slip', 'dependent', 'Moment' );

%% Optimization Options
Fit.Opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
Fit.Opts.Display = 'Off';

Fit.Opts.MaxFunEvals = 10000;
Fit.Opts.MaxIter = 10000;

Fit.Opts.StartPoint = [0.25 1                  0.005 0.1 0.25  0  1];
Fit.Opts.Lower =      [   0 0 -max(abs(Data.Moment))   0   -5 -5  0];
Fit.Opts.Upper =      [   2 5  max(abs(Data.Moment))   1    1  5 15];

%% Fit model to data
[Fit.Result, Fit.GoF] = fit( Data.Slip, Data.Moment, Fit.Type, Fit.Opts );

%% Allocate Results
Nominal.Bt = Fit.Result.Bt;
Nominal.Ct = Fit.Result.Ct;
Nominal.Dt = Fit.Result.Dt;
Nominal.Et = Fit.Result.Et;
Nominal.Ht = Fit.Result.Ht;

Nominal.qbz10 = Fit.Result.qbz10;
Nominal.Cr = 1;
Nominal.Dr = Fit.Result.Dr;
Nominal.Hf = Fyo.Hy + Fyo.Vy/Fyo.Kya;

Nominal.Residual = ( Data.Moment - Fit.Result( Data.Slip ) )';

%{
figure( 'Name', 'Aligning Moment Test', 'NumberTitle', 'off' )
plot( Data.Slip, Data.Moment, 'k.', Data.Slip, Fit.Result( Data.Slip ), 'g-.', Data.Slip, Nominal.Residual, 'r.' )
pause(0.25)
delete( gcf )
%}

%% Local Functions
    function Fyo = FyoEvaluation( Mesh, Tire )
        % Note Function is evaluated for null camber
        Fyo.Cy = Tire.Pacejka.p.C.y(1);
        
        Fyo.Dy = (Tire.Pacejka.p.D.y(1) + Tire.Pacejka.p.D.y(2).*Mesh.dFz) .* ...
            (1 + Tire.Pacejka.p.P.y(3).*Mesh.dPi + Tire.Pacejka.p.P.y(4).*Mesh.dPi.^2) .* ...
            (1 - Tire.Pacejka.p.D.y(3).*0.^2).*Mesh.Load;
        
        Fyo.Ey(1) = ( Tire.Pacejka.p.E.y(1) + Tire.Pacejka.p.E.y(2).*Mesh.dFz ) .* ...
            ( 1 + Tire.Pacejka.p.E.y(5).*0.^2 + ...
            ( Tire.Pacejka.p.E.y(3) + Tire.Pacejka.p.E.y(4).*0 ) );
        
        Fyo.Ey(2) = ( Tire.Pacejka.p.E.y(1) + Tire.Pacejka.p.E.y(2).*Mesh.dFz ) .* ...
            ( 1 + Tire.Pacejka.p.E.y(5).*0.^2 - ...
            ( Tire.Pacejka.p.E.y(3) + Tire.Pacejka.p.E.y(4).*0 ) );
        
        Fyo.Kya = Tire.Pacejka.p.K.y(1) .* Tire.Pacejka.Fzo .* ( 1 + Tire.Pacejka.p.P.y(1).*Mesh.dPi ) .* ...
            ( 1 - Tire.Pacejka.p.K.y(3).*abs(0) ) .* sin( Tire.Pacejka.p.K.y(4) .* ...
            atan( (Mesh.Load./Tire.Pacejka.Fzo) ./ ...
            ( ( Tire.Pacejka.p.K.y(2) + Tire.Pacejka.p.K.y(5).*0.^2 ) .* ( 1 + Tire.Pacejka.p.P.y(2).*Mesh.dPi ) ) ) );
        
        Fyo.Kyg0 = Mesh.Load.*(Tire.Pacejka.p.K.y(6) + Tire.Pacejka.p.K.y(7).*Mesh.dFz ) .* (1 + Tire.Pacejka.p.P.y(5).*Mesh.dPi );
        
        Fyo.By = Fyo.Kya ./ ( Fyo.Cy.*Fyo.Dy );
        
        Fyo.Vyg = Mesh.Load.*(Tire.Pacejka.p.V.y(3) + Tire.Pacejka.p.V.y(4).*Mesh.dFz ).*0;
        
        Fyo.Vy = Mesh.Load.*(Tire.Pacejka.p.V.y(1) + Tire.Pacejka.p.V.y(2).*Mesh.dFz ) + Fyo.Vyg;
        
        Fyo.Hy = (Tire.Pacejka.p.H.y(1) + Tire.Pacejka.p.H.y(2).*Mesh.dFz ) .* ...
            (Fyo.Kyg0.*0 - Fyo.Vyg ) ./ Fyo.Kya;
    end
end
