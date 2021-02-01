 function [Nominal] = PureAligningNominal( Mesh, Raw, Tire )
%% Evaluating Pure Lateral Force Fit (Fyo)
[Fyo, By, Cy, Kya, Hy, Vy] = FyoEval( Tire );

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

% Coefficients           Bt  Ct                     Dr  Dt   Et Ht qbz10
Fit.Opts.StartPoint = [0.25   1                  0.005 0.1 0.25  0   1];
Fit.Opts.Lower =      [   0   1 -max(abs(Data.Moment))   0  -10 -5   0];
Fit.Opts.Upper =      [  50 1.1  max(abs(Data.Moment)) Inf    1  5  50];

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

figure( 'Name', 'Aligning Moment Test', 'NumberTitle', 'off' )
plot( Data.Slip, Data.Moment, 'k.', Data.Slip, Fit.Result( Data.Slip ), 'g-.', Data.Slip, Nominal.Residual, 'r.' )
pause(0.25)
delete( gcf )

%% Local Functions
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
            ( Tire.Pacejka.p.E.y(3) + Tire.Pacejka.p.E.y(4).*Gam ).*sign(deg2rad(Slip) + Hy) );

        Fyo = @(Pi, Fz, Gam, Slip) Dy(Pi, Fz, Gam) .* ...
            sin( Cy .* atan( (1-Ey(Fz, Gam, Slip, Hy(Pi, Fz, Gam) )) .* ...
            By(Pi, Fz, Gam).*(deg2rad(Slip) + Hy(Pi, Fz, Gam) ) + ...
            Ey(Fz, Gam, Slip, Hy(Pi, Fz, Gam) ).*atan( ...
            By(Pi, Fz, Gam).*(deg2rad(Slip) + Hy(Pi, Fz, Gam) ) ) ) ) + Vy(Fz, Gam);
    end
end
