function [Fx, Fy, Mz, Mx, My] = ContactPatchLoads( Tire, ...
    SlipAngle, SlipRatio, ...
    NormalLoad, Pressure, Inclination, Velocity, ...
    Idx, Model )

%% ContactPatchLoads - Tire Load Evaluation
% This script evaluates a Pacejka contact patch tire model using the
% formatted by TireModelFittingMain(). 
% 
% Inputs:
%   Tire          - (struct)      Pacejka Parameters
%   SlipAngle     - (n,1 numeric) Slip Angle          {alpha} [deg]
%   SlipRatio     - (n,1 numeric) Slip Ratio          {kappa} [ ]
%   NormalLoad    - (n,1 numeric) Normal Load         {F_z}   [N]
%   Pressure      - (n,1 numeric) Inflation Pressure  {P_i}   [kPa]
%   Inclination   - (n,1 numeric) Inclination Angle   {gamma} [deg]
%   Velocity      - (n,1 numeric) Center Velocity     {v_c}   [m/s]
%   Idx           - (n,1 numeric) Corner Idx          {1-4}   [ ]
%   Model         - (struct)      Fidelity Choices 
%       .Pure     - (char)        Pure Slip Model     {'Linear', 'Pacejka'}
%       .Combined - (char)        Combined Slip Model {'Pure', 'MNC', 'Pacejka'}
%
% Outputs:
%   Fx - (n,1 numeric) Longitudinal Force {F_x} [N]
%   Fy - (n,1 numeric) Lateral Force      {F_y} [N]
%   Mz - (n,1 numeric) Aligning Moment    {M_z} [N-m]
%   Mx - (n,1 numeric) Overturning Moment {M_x} [N-m]
%   My - (n,1 numeric) Rolling Resistance {M_y} [N-m]
%
% Notes:
%   Pacejka combined slip model is not current implemented (4/2/21) and
%   when the model does become available the model compatibility checks
%   should be reformatted.
%
% Author(s): 
% Blake Christierson (bechristierson@ucdavis.edu) [Sep 2018 - Jun 2021] 
% Carlos Lopez       (calopez@ucdavis.edu       ) [Jan 2019 -         ]
% 
% Last Updated: 2-Apr-2021

%%% Model Compatibility Check
if ~ismember(Model.Pure, {'Linear', 'Pacejka'})
    error('Please choose either ''Linear'' or ''Pacejka'' for the pure slip model')
elseif ~ismember(Model.Combined, {'Pure', 'MNC'})
    if strcmp(Model.Combined, 'Pacejka')
        error('The Pacejka combined slip model is not implemented (4/2/21)')
    else
        error('Please choose either ''Pure'' or ''MNC'' for the combined slip model')
    end
end

%%% Slip Angle Conversion
SlipAngle = deg2rad(SlipAngle);

%%% Left-Right Side Orientation Manipulation
% Ensures model conicity & ply-steer are symmetric to prevent drift
SlipAngle = SlipAngle .* (-1).^(mod(Idx+1,2));
Inclination   = Inclination   .* (-1).^(mod(Idx+1,2));

%%% Nondimensional Load & Pressure
dFz = (NormalLoad - Tire.Pacejka.Fzo) ./ Tire.Pacejka.Fzo;
dPi = (Pressure - Tire.Pacejka.Pio) ./ Tire.Pacejka.Pio;

%%% Evaluate Pure Slip Forces & Aligning Moment
[Fx0, Kxk, Kappa0] = EvaluateFx0( Tire, SlipAngle, SlipRatio, ...
    NormalLoad, dFz, Pressure, dPi, Inclination, Velocity, Idx, Model ); 

[Fy0, By, Cy, Kya, Hy, Vy, Alpha0] = EvaluateFy0( Tire, SlipAngle, SlipRatio, ...
    NormalLoad, dFz, Pressure, dPi, Inclination, Velocity, Idx, Model ); 
    
[Fyo, ~, ~, ~, ~, ~, Alphao] = EvaluateFy0( Tire, SlipAngle, SlipRatio, ...
    NormalLoad, dFz, Pressure, dPi, 0  , Velocity, Idx, Model );

[Mz0, Bt, Ct, Dt, Et, Ht, Br, Cr, Dr, Hf] = EvaluateMz0( Tire, SlipAngle, SlipRatio, ...
    NormalLoad, dFz, Pressure, dPi, Inclination, Velocity, Idx, Model, ...
    Fyo, By, Cy, Kya, Hy, Vy );

%%% Combined Slip Force Models 
switch Model.Combined
    case 'Pure' % Pure Slip Only Model
        Fx  = Fx0;
        Fy  = Fy0;   

    case 'MNC' % Modified-Nicolas-Comstock Combined Slip Model
        [Fx, Fy ]  = MNCCombinedSlipForces( Tire, ...
            SlipAngle, SlipRatio, Fx0, Fy0, Kxk, Kya, Kappa0, Alpha0 );
        [~ , Fyp]  = MNCCombinedSlipForces( Tire, ...
            SlipAngle, SlipRatio, Fx0, Fyo, Kxk, Kya, Kappa0, Alphao );

    case 'MF6_1' % Combined Slip Pacejka MF6.1 Model
        % Not Implemented
end

%%% Evaluate Combined Slip Aligning Moment (Mz)
if strcmpi( Model.Combined, 'Pure' )
    Mz = Mz0;
else
    Mz = EvaluateMz( Tire, SlipAngle, SlipRatio, ...
        NormalLoad, dFz, Pressure, dPi, Inclination, Velocity, Idx, Model, ...
        Fx, Fy, Fyp, Kxk, Kya, Bt, Ct, Dt, Et, Ht, Br, Cr, Dr, Hf );
end

%%% Evaluate Moment & Rolling Resistance
Mx = EvaluateMx( Tire, SlipAngle, SlipRatio, ...
    NormalLoad, dFz, Pressure, dPi, Inclination, Velocity, Idx, Model, Fy );

My = EvaluateMy( Tire, SlipAngle, SlipRatio, ...
    NormalLoad, dFz, Pressure, dPi, Inclination, Velocity, Idx, Model, Fx );

%%% Pure Slip Longitudinal Force Evaluation (Fx0)
function [Fx0, Kxk, Kappa0] = ...
        EvaluateFx0( Tire, ~, Kappa, Fz, dFz, ~, dPi, Inc, ~, ~, Model )
    % Evaluate P6 Pacejka
    Cx = Tire.Pacejka.p.C.x(1) .* Tire.Pacejka.L.C.x;

    Dx = (Tire.Pacejka.p.D.x(1) + Tire.Pacejka.p.D.x(2).*dFz) .* ...
         (1 + Tire.Pacejka.p.P.x(3).*dPi + Tire.Pacejka.p.P.x(4).*dPi.^2) .* ...
         (1 - Tire.Pacejka.p.D.x(3).*Inc.^2).*Fz .* Tire.Pacejka.L.mu.x;

    Ex = ( Tire.Pacejka.p.E.x(1) + Tire.Pacejka.p.E.x(2).*dFz ...
        + Tire.Pacejka.p.E.x(3).*dFz.^2 ) .* ( 1 - Tire.Pacejka.p.E.x(4).*sign(Kappa) ) .* ...
        Tire.Pacejka.L.E.x;

    Kxk = Fz.*(Tire.Pacejka.p.K.x(1) + Tire.Pacejka.p.K.x(2).*dFz ) .* ...
        exp( Tire.Pacejka.p.K.x(3) .* dFz ) .* ...
        (1 + Tire.Pacejka.p.P.x(1).*dPi + Tire.Pacejka.p.P.x(2).*dPi.^2) .* ...
        Tire.Pacejka.L.K.x.k;

    Bx = Kxk ./ ( Cx.*Dx );

    Vx = Fz.*(Tire.Pacejka.p.V.x(1) + Tire.Pacejka.p.V.x(2).*dFz ) .* ...
        Tire.Pacejka.L.V.x;

    Hx = (Tire.Pacejka.p.H.x(1) + Tire.Pacejka.p.H.x(2).*dFz) .* ...
        Tire.Pacejka.L.H.x;

    % Evaluate Longitudinal Force​
    switch Model.Pure
        case 'Linear'
            Fx0 = Kxk .* Kappa;
        case 'Pacejka'
            Fx0 = Dx .* sin( Cx .* atan( (1-Ex) .* Bx.*(Kappa + Hx) + ...
                Ex.*atan( Bx.*(Kappa + Hx) ) ) ) + Vx;
    end

    % Evaluate Null Slip Ratio
    if strcmpi( Model.Combined, 'MNC' )
        Kappa0 = ((Fx0 .* Hx)./Kxk) + Vx;
%                  Opts = optimoptions( 'fsolve', 'Display', 'off', 'Algorithm', 'Levenberg-Marquardt' );
%                 Kappa0 = fsolve( @(k) Dx .* sin( Cx .* atan( (1-Ex) .* Bx.*(k + Hx) + ...
%                      Ex.*atan( Bx.*(k + Hx) ) ) ) + Vx, zeros(size(Dx)), Opts );
    else
        Kappa0 = [];
    end
end

%%% Pure Slip Lateral Force Evaluation (Fy0)
function [Fy0, By, Cy, Kya, Hy, Vy, Alpha0] = ...
        EvaluateFy0( Tire, Alpha, ~, Fz, dFz, ~, dPi, Inc, ~, i, Model )            
    % Evaluate P6 Pacejka
    Cy = Tire.Pacejka.p.C.y(1);

    Dy = (Tire.Pacejka.p.D.y(1) + Tire.Pacejka.p.D.y(2).*dFz) .* ...
         (1 + Tire.Pacejka.p.P.y(3).*dPi + Tire.Pacejka.p.P.y(4).*dPi.^2) .* ...
         (1 - Tire.Pacejka.p.D.y(3).*Inc.^2).*Fz;

    Kya = Tire.Pacejka.p.K.y(1) .* Tire.Pacejka.Fzo .* ( 1 + Tire.Pacejka.p.P.y(1).*dPi ) .* ...
        ( 1 - Tire.Pacejka.p.K.y(3).*abs(Inc) ) .* sin( Tire.Pacejka.p.K.y(4) .* ...
        atan( (Fz./Tire.Pacejka.Fzo) ./ ...
        ( ( Tire.Pacejka.p.K.y(2) + Tire.Pacejka.p.K.y(5).*Inc.^2 ) .* ...
        ( 1 + Tire.Pacejka.p.P.y(2).*dPi ) ) ) );

    Kyg0 = Fz.*(Tire.Pacejka.p.K.y(6) + Tire.Pacejka.p.K.y(7).*dFz ) .* ...
        (1 + Tire.Pacejka.p.P.y(5).*dPi );

    By = Kya ./ ( Cy.*Dy );

    Vyg = Fz.*(Tire.Pacejka.p.V.y(3) + Tire.Pacejka.p.V.y(4).*dFz ).*Inc;

    Vy = Fz.*(Tire.Pacejka.p.V.y(1) + Tire.Pacejka.p.V.y(2).*dFz ) + Vyg;

    Hy = (Tire.Pacejka.p.H.y(1) + Tire.Pacejka.p.H.y(2).*dFz ) .* (Kyg0.*Inc - Vyg ) ./ Kya;

    Ey = ( Tire.Pacejka.p.E.y(1) + Tire.Pacejka.p.E.y(2).*dFz ) .* ...
        ( 1 + Tire.Pacejka.p.E.y(5).*Inc.^2 - ...
        ( Tire.Pacejka.p.E.y(3) + Tire.Pacejka.p.E.y(4).*Inc ).*sign(Alpha + Hy) );

    % Evaluate Lateral Force
    switch Model.Pure
        case 'Linear'
            Fy0 = Kya .* Alpha;
        case 'Pacejka'
            Fy0 = (Dy .* sin( Cy .* atan( (1-Ey) .* By.*( Alpha + Hy ) + ...
                Ey .* atan( By.*( Alpha + Hy ) ) ) ) + Vy) .* (-1).^(mod(i+1,2));
    end

    % Evaluate Null Slip Angle
    if strcmpi( Model.Combined, 'MNC' )
        Alpha0 = ((Fy0 .* Hy)./Kya) + Vy;
%                  Opts = optimoptions( 'fsolve', 'Display', 'off', 'Algorithm', 'Levenberg-Marquardt' );
%                    Alpha0 = fsolve( @(a) (Dy .* sin( Cy .* atan( (1-Ey) .* By.*(a +Hy) + ...
%                           Ey.*atan( By.*(a + Hy) ) ) ) + Vy) .* (-1).^(mod(i+1,2)), zeros(size(Dy)), Opts );

    else
        Alpha0 = [];
    end
end

%%% Evaluate Pure Slip Aligning Moment (Mz0)
function [Mz0, Bt, Ct, Dt, Et, Ht, Br, Cr, Dr, Hf] = ...
        EvaluateMz0( Tire, Alpha, ~, Fz, dFz, ~, dPi, Inc, ~, i, ~, ...
            Fyo, By, Cy, Kya, Hy, Vy )
    % Evaluate P6 Pacejka
    Bt = (Tire.Pacejka.q.B.z(1) + Tire.Pacejka.q.B.z(2).*dFz + ...
        Tire.Pacejka.q.B.z(3).*dFz.^2) .* ...
        (1 + Tire.Pacejka.q.B.z(5).*abs(Inc) + ...
        Tire.Pacejka.q.B.z(6).*Inc.^2);

    Ct = Tire.Pacejka.q.C.z(1);

    Dt = (Tire.Pacejka.Ro ./ Tire.Pacejka.Fzo) .* ...
        ( Tire.Pacejka.q.D.z(5) .* (Tire.Pacejka.Fzo ./ Tire.Pacejka.Ro) + ...
        Tire.Pacejka.q.D.z(1).*Fz + Tire.Pacejka.q.D.z(2).*Fz.*dFz ) .* ...
        ( 1 - Tire.Pacejka.p.P.z(1).*dPi ) .* ...
        ( 1 + Tire.Pacejka.q.D.z(3).*abs(Inc) + ...
        Tire.Pacejka.q.D.z(4).*Inc.^2 );

    Ht = Tire.Pacejka.q.H.z(1) + ...
        Tire.Pacejka.q.H.z(2).*dFz + (Tire.Pacejka.q.H.z(3) + ...
        Tire.Pacejka.q.H.z(4).*dFz).*Inc;

    Et = (Tire.Pacejka.q.E.z(1) + Tire.Pacejka.q.E.z(2).*dFz + ...
        Tire.Pacejka.q.E.z(3).*dFz.^2) .* ...
        (1 + (Tire.Pacejka.q.E.z(4) + Tire.Pacejka.q.E.z(5).*Inc).*(2/pi).* ...
        atan( Bt .* Ct .* (Alpha+Ht) ) );

    Br = Tire.Pacejka.q.B.z(10) .* By .* Cy;

    Cr = 1;

    Dr = Tire.Pacejka.Ro .* Fz .* ( ( Tire.Pacejka.q.D.z(6) + ...
        Tire.Pacejka.q.D.z(7).*dFz ) + ...
        ( ( Tire.Pacejka.q.D.z(8) + Tire.Pacejka.q.D.z(9).*dFz ) .* ...
        ( 1 + Tire.Pacejka.p.P.z(2).*dPi ) + ...
        ( Tire.Pacejka.q.D.z(10) + Tire.Pacejka.q.D.z(11).*dFz ) .* ...
        abs(Inc) ).*Inc ) .* cos( Alpha );

    Hf = Hy + Vy ./ Kya;

    % Evaluate Trail Function
    t0 = Dt.*cos( Ct.*atan( (1-Et).*(Bt.*(Alpha + Ht)) + ...
        Et.*atan( Bt.*(Alpha + Ht) ) ) ) .* cos( Alpha ) ;

    % Evaluate Residual Aligning Moment
    Mzro = Dr .* cos( Cr.*atan( Br.*(Alpha + Hf) ) ).* cos( Alpha );

    % Evaluate Pure Slip Aligning Moment (Mz0)
    Mz0 = (-t0 .* Fyo .* (-1).^(mod(i+1,2)) + Mzro) .* (-1).^(mod(i+1,2));
end

%%% Modified-Nicolas Comstock Combined Slip Forces
function [Fx, Fy] = MNCCombinedSlipForces( ~, Alpha, Kappa, Fx0, Fy0, Kxk, Kya, Kappa0, Alpha0 )
    Fx = abs(Fx0 .* Fy0 ./ sqrt( (Kappa - Kappa0).^2 .* Fy0.^2 + Fx0.^2 .* tan(Alpha - Alpha0).^2 ) .* ...
        sqrt( (Kappa - Kappa0).^2 .* Kya.^2 + (1 - abs(Kappa - Kappa0)).^2 .* cos(Alpha - Alpha0).^2 .* Fx0.^2 ) ./ Kya) .* sign(Fx0);

    Fy = abs( Fx0 .* Fy0 ./ sqrt( (Kappa - Kappa0).^2 .* Fy0.^2 + Fx0.^2 .* tan(Alpha - Alpha0).^2 ) .* ...
        sqrt( (1 - abs(Kappa - Kappa0)).^2 .* cos(Alpha - Alpha0).^2 .* Fy0.^2 + sin(Alpha - Alpha0).^2 .* Kxk.^2 ) ./ ...
        ( Kxk .* cos(Alpha - Alpha0)) ) .* sign(Fy0);
end

%%% Evaluate Combined Slip Aligning Moment (Mz)
function Mz = EvaluateMz( Tire, Alpha, Kappa, ~, dFz, ~, ~, Inc, ~, i, ~, ...
    Fx, Fy, Fyp, Kxk, Kya, Bt, Ct, Dt, Et, Ht, Br, Cr, Dr, Hf )
    % Equivalent Slips
    Alphat = (Alpha + (-1).^(mod(i+1,2)).*Ht );
    Alphar = (Alpha + (-1).^(mod(i+1,2)).*Hf );

    Alphateq = sqrt( Alphat.^2 + (Kxk./Kya).^2 .* Kappa.^2 ) .* sign( Alphat );
    Alphareq = sqrt( Alphar.^2 + (Kxk./Kya).^2 .* Kappa.^2 ) .* sign( Alphar );

    % Evaluate Trail Function
    t = Dt.*cos( Ct.*atan( (1-Et).*(Bt.*Alphateq) + ...
        Et.*atan(Bt.*Alphateq) ) ) .* cos( Alpha );

    % Evaluate Scrub Function
    s = Tire.Pacejka.Ro .* (Tire.Pacejka.s.s.z(1) + ...
        Tire.Pacejka.s.s.z(2) .* (Fy./Tire.Pacejka.Fzo) + ...
        (Tire.Pacejka.s.s.z(3) + Tire.Pacejka.s.s.z(4) .* dFz) .* Inc );

    % Evaluate Residual Aligning Moment
    Mzr = Dr .* cos( Cr.*atan( Br.*(Alphareq) ) ) .* cos( Alpha );

    % Evaluate Pure Slip Aligning Moment (Mz0)
    Mz = -t .* Fyp + Mzr + s .* Fx;
end

%%% Evaluate Overturning Moment
function Mx = EvaluateMx( Tire, ~, ~, Fz, ~, ~, dPi, Inc, ~, ~, ~, Fy )
    Mx = 0;
    %{
    Mx = Tire.Pacejka.Ro .* Fz .* ( Tire.Pacejka.q.s.x(1) - ...
        Tire.Pacejka.q.s.x(2) .* Inc .* ( 1 + Tire.Pacejka.p.P.Mx(1).*dPi ) + ...
        Tire.Pacejka.q.s.x(3) .* Fy ./ Tire.Pacejka.Fzo + ...
        Tire.Pacejka.q.s.x(4) .* cos( Tire.Pacejka.q.s.x(5) .* ...
        atan( (Tire.Pacejka.q.S.x(6) .* Fz./Tire.Pacejka.Fzo) ).^2 ) .* ...
        sin( Tire.Pacejka.q.s.x(7) .* Inc + Tire.Pacejka.q.s.x(8) .* ...
        atan( Tire.Pacejka.q.s.x(9) .* Fy./Tire.Pacejka.Fzo ) ) + ...
        Tire.Pacejka.q.s.x(10) .* atan( Tire.Pacejka.q.s.x(11) .* Fz./Tire.Pacejka.Fzo) .* Inc );
    %}
end

%%% Evaluate Rolling Resistance
function My = EvaluateMy( Tire, Alpha, Kappa, Fz, dFz, Pi, dPi, Inc, V, ~, Model, Fx )
    My = 0;
end

end