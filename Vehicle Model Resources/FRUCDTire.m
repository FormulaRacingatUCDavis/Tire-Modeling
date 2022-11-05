%% FRUCD Wheel & Tire Class Definition
% Defines wheel & tire parameters for use in vehicle simulation and
% provides methods for evaluating certain properties. The tire operates in
% the SAE Z-up coordinate system.
%
% Main Methods:
% - FRUCDTire         : Instance Generation
% - ContactPatchLoads : Evaluates Model for F_x, F_y, M_z, M_x, and M_y
% - RadialDisplacement: Evaluates Loaded & Effective Radius (r_l & r_e)
%
% Variables:
% - Slip Angle        : alpha [deg]
% - Slip Ratio        : kappa [   ]
% - Inflation Pressure: P_i   [kPa]
% - Inclination Angle : gamma [deg]
% - Tire Velocity     : v     [m/s]
% - Normal Load       : F_z   [N  ]
% - Longitudinal Force: F_x   [N  ]
% - Lateral Force     : F_y   [N  ]
% - Aligning Moment   : M_z   [Nm ]
% - Overturning Moment: M_x   [Nm ]
% - Rolling Resistance: M_y   [Nm ]
% - Effective Radius  : r_e   [mm ]
% - Loaded Radius     : r_l   [mm ]
%
% Author(s):
% Blake Christierson (bechristierson@ucdavis.edu)
%
% Last Updated: 14-Feb-2021

classdef FRUCDTire
    properties
        Name
        Date
        Time
        Source
        Notes
        Pacejka
        Radius
        Mass
        Thermal
    end
    methods
        function obj = FRUCDTire(Name, Source, Notes)
            if nargin > 0
                obj.Name    = Name;                      % Tire Name
                obj.Date    = datestr(now,'dd-mm-yyyy'); % Generation Date
                obj.Time    = datestr(now,'HH:MM PM'  ); % Generation Time
                obj.Source  = Source;                    % Data Sources     
                obj.Notes   = Notes;                     % Notes About Model
                
                obj.Pacejka = PacejkaParameters();       % Pacejka Parameters
                obj.Radius  = RadiusParameters();        % Radius Parameters
                obj.Mass    = MassParameters();          % Mass Parameters
                obj.Thermal = ThermalParameters();       % Thermal Parameters
            end
            
            function Pacejka = PacejkaParameters()
                %%% Nominal Conditions
                Pacejka.Fzo = 3000;      % Nominal Pacejka Load     [N] 
                Pacejka.Pio = 70;        % Nominal Pacejka Pressure [kPa]
                Pacejka.Ro = 8*2.54/100; % Nominal Pacejka Radius   [m]

                % Scaling Factors (Lambdas)
                Pacejka.L.F.zo = 1; 

                Pacejka.L.mu.x = 1;
                Pacejka.L.mu.y = 1;
                Pacejka.L.mu.V = 0;

                Pacejka.L.K.x.k = 1;
                Pacejka.L.K.x.a = 1;

                Pacejka.L.K.y.a = 1;
                Pacejka.L.K.y.k = 1;

                Pacejka.L.K.y.g = 1;
                Pacejka.L.K.z.g = 1;

                Pacejka.L.C.x = 1;
                Pacejka.L.C.y = 1;
                Pacejka.L.C.z = 1;

                Pacejka.L.E.x = 1;
                Pacejka.L.E.y = 1;

                Pacejka.L.H.x = 1;
                Pacejka.L.H.y = 1;

                Pacejka.L.V.x = 1;
                Pacejka.L.V.y = 1;
                Pacejka.L.V.yk = 1;

                Pacejka.L.t = 1;

                Pacejka.L.M.x = 1;
                Pacejka.L.M.y = 1;
                Pacejka.L.M.r = 1;

                % Pacejka Spin Factors (Zetas)
                Pacejka.Z = ones(8,1);

                % p-Factors (Pure Slip Force Coefficients)
                Pacejka.p.C.x = 0;
                Pacejka.p.D.x = zeros(3,1);
                Pacejka.p.E.x = zeros(4,1);
                Pacejka.p.K.x = zeros(3,1);
                Pacejka.p.H.x = zeros(2,1);
                Pacejka.p.V.x = zeros(2,1);
                Pacejka.p.P.x = zeros(4,1);

                Pacejka.p.C.y = 0;
                Pacejka.p.D.y = zeros(3,1);
                Pacejka.p.E.y = zeros(5,1);
                Pacejka.p.K.y = zeros(7,1);
                Pacejka.p.H.y = zeros(2,1);
                Pacejka.p.V.y = zeros(4,1);
                Pacejka.p.P.y = zeros(5,1);

                % q-Factors (Moment Coefficients)
                Pacejka.q.B.z = zeros(10,1);
                Pacejka.q.C.z = 0;
                Pacejka.q.D.z = zeros(11,1);
                Pacejka.q.E.z = zeros(5,1);
                Pacejka.q.H.z = zeros(4,1);
                Pacejka.p.P.z = zeros(2,1);

                Pacejka.q.s.x = zeros(10,1);
                Pacejka.q.s.y = zeros(8,1);

                % r-Factors (Combined Slip Force Coefficients)
                Pacejka.r.B.x = zeros(3,1);
                Pacejka.r.C.x = 0;
                Pacejka.r.E.x = zeros(2,1);
                Pacejka.r.H.x = 0;

                Pacejka.r.B.y = zeros(4,1);
                Pacejka.r.C.y = 0;
                Pacejka.r.E.y = zeros(2,1);
                Pacejka.r.H.y = zeros(2,1);
                Pacejka.r.V.y = zeros(6,1);

                % s-Factors (Combined Slip Aligning Moment Coefficients)
                Pacejka.s.s.z = zeros(4,1);
            end
            
            function Radius = RadiusParameters()
                Radius.Effective = [];
                Radius.Loaded    = [];
            end
            
            function Mass = MassParameters()
                Mass.m  = 0;            % Wheel Package Mass [kg]
                Mass.I  = zeros(3,3);   % Wheel Pacakge Inertia Tensor [kg-m^2]
                
                Mass.Is = 0;            % Spin Inertia (Spindle, Rotor, Rim, ...) [kg-m^2]
            end
            
            function Thermal = ThermalParameters()
                Thermal = [];
            end
        end    
        
        %%% Evaluate Tire Model for Fx, Fy, Mz, Mx, My
        function [Fx, Fy, Mz, Mx, My] = ContactPatchLoads( obj, Alpha, Kappa, Fz, Pi, Inc, V, i, Fidelity )   
            %%% Slip Angle Conversion
            Alpha = deg2rad(Alpha);
            
            %%% Side Manipulation
            Alpha = Alpha .* (-1).^(mod(i+1,2));
            Inc   = Inc   .* (-1).^(mod(i+1,2));
            
            %%% Nondimensional Load & Pressure
            dFz = (Fz - obj.Pacejka.Fzo) ./ obj.Pacejka.Fzo;
            dPi = (Pi - obj.Pacejka.Pio) ./ obj.Pacejka.Pio;
            
            %%% Evaluate Pure Slip Forces & Aligning Moment
            [Fx0, Kxk, Kappa0] = ...
                EvaluateFx0( obj, Alpha, Kappa, Fz, dFz, Pi, dPi, Inc, V, i, Fidelity ); 
            
            [Fy0, By, Cy, Kya, Hy, Vy, Alpha0] = ...
                EvaluateFy0( obj, Alpha, Kappa, Fz, dFz, Pi, dPi, Inc, V, i, Fidelity ); 
            
            [Fyo, ~, ~, ~, ~, ~, Alphao] = ...
                EvaluateFy0( obj, Alpha, Kappa, Fz, dFz, Pi, dPi, 0  , V, i, Fidelity );
            
            [Mz0, Bt, Ct, Dt, Et, Ht, Br, Cr, Dr, Hf] = ...
                EvaluateMz0( obj, Alpha, Kappa, Fz, dFz, Pi, dPi, Inc, V, i, Fidelity, ...
                    Fyo, By, Cy, Kya, Hy, Vy );
            
            %%% Combined Slip Force Models 
            switch Fidelity.Combined
                case 'Pure' % Pure Slip Only Model
                    Fx  = Fx0;
                    Fy  = Fy0;   
                    
                case 'MNC' % Modified-Nicolas-Comstock Combined Slip Model
                    [Fx, Fy ] = MNCCombinedSlipForces( obj, ...
                        Alpha, Kappa, Fx0, Fy0, Kxk, Kya, Kappa0, Alpha0 );
                    [~ , Fyp]  = MNCCombinedSlipForces( obj, ...
                        Alpha, Kappa, Fx0, Fyo, Kxk, Kya, Kappa0, Alphao );
                         
                case 'MF6_1' % Combined Slip Pacejka MF6.1 Model
                    % Not Implemented
            end
            
            %%% Evaluate Combined Slip Aligning Moment (Mz)
            if strcmpi( Fidelity.Combined, 'Pure' )
                Mz = Mz0;
            else
                Mz = EvaluateMz( obj, Alpha, Kappa, Fz, dFz, Pi, dPi, Inc, V, i, Fidelity, ...
                    Fx, Fy, Fyp, Kxk, Kya, Bt, Ct, Dt, Et, Ht, Br, Cr, Dr, Hf );
            end
            
            %%% Evaluate Moment & Rolling Resistance
            Mx = EvaluateMx( obj, Alpha, Kappa, Fz, dFz, Pi, dPi, Inc, V, i, Fidelity, Fy );
            
            My = EvaluateMy( obj, Alpha, Kappa, Fz, dFz, Pi, dPi, Inc, V, i, Fidelity, Fx );
            
            %%% Pure Slip Longitudinal Force Evaluation (Fx0)
            function [Fx0, Kxk, Kappa0] = ...
                    EvaluateFx0( obj, ~, Kappa, Fz, dFz, ~, dPi, Inc, ~, ~, Fidelity )
                % Evaluate P6 Pacejka
                Cx = obj.Pacejka.p.C.x(1) .* obj.Pacejka.L.C.x;

                Dx = (obj.Pacejka.p.D.x(1) + obj.Pacejka.p.D.x(2).*dFz) .* ...
                     (1 + obj.Pacejka.p.P.x(3).*dPi + obj.Pacejka.p.P.x(4).*dPi.^2) .* ...
                     (1 - obj.Pacejka.p.D.x(3).*Inc.^2).*Fz .* obj.Pacejka.L.mu.x;

                Ex = ( obj.Pacejka.p.E.x(1) + obj.Pacejka.p.E.x(2).*dFz ...
                    + obj.Pacejka.p.E.x(3).*dFz.^2 ) .* ( 1 - obj.Pacejka.p.E.x(4).*sign(Kappa) ) .* ...
                    obj.Pacejka.L.E.x;

                Kxk = Fz.*(obj.Pacejka.p.K.x(1) + obj.Pacejka.p.K.x(2).*dFz ) .* ...
                    exp( obj.Pacejka.p.K.x(3) .* dFz ) .* ...
                    (1 + obj.Pacejka.p.P.x(1).*dPi + obj.Pacejka.p.P.x(2).*dPi.^2) .* ...
                    obj.Pacejka.L.K.x.k;

                Bx = Kxk ./ ( Cx.*Dx );

                Vx = Fz.*(obj.Pacejka.p.V.x(1) + obj.Pacejka.p.V.x(2).*dFz ) .* ...
                    obj.Pacejka.L.V.x;

                Hx = (obj.Pacejka.p.H.x(1) + obj.Pacejka.p.H.x(2).*dFz) .* ...
                    obj.Pacejka.L.H.x;

                % Evaluate Longitudinal Force​
                switch Fidelity.Pure
                    case 'Linear'
                        Fx0 = Kxk .* Kappa;
                    case 'Pacejka'
                        Fx0 = Dx .* sin( Cx .* atan( (1-Ex) .* Bx.*(Kappa + Hx) + ...
                            Ex.*atan( Bx.*(Kappa + Hx) ) ) ) + Vx;
                end

                % Evaluate Null Slip Ratio
                if strcmpi( Fidelity.Combined, 'MNC' )
                    Kappa0 = -Vx ./ Kxk - Hx;
 %                   Opts = optimoptions( 'fsolve', 'Display', 'off', 'Algorithm', 'Levenberg-Marquardt' );
  %                  Kappa0 = fsolve( @(k) Dx .* sin( Cx .* atan( (1-Ex) .* Bx.*(k + Hx) + ...
   %                     Ex.*atan( Bx.*(k + Hx) ) ) ) + Vx, zeros(size(Dx)), Opts );
                else
                    Kappa0 = [];
                end
            end

            %%% Pure Slip Lateral Force Evaluation (Fy0)
            function [Fy0, By, Cy, Kya, Hy, Vy, Alpha0] = ...
                    EvaluateFy0( obj, Alpha, ~, Fz, dFz, ~, dPi, Inc, ~, i, Fidelity )            
                % Evaluate P6 Pacejka
                Cy = obj.Pacejka.p.C.y(1);

                Dy = (obj.Pacejka.p.D.y(1) + obj.Pacejka.p.D.y(2).*dFz) .* ...
                     (1 + obj.Pacejka.p.P.y(3).*dPi + obj.Pacejka.p.P.y(4).*dPi.^2) .* ...
                     (1 - obj.Pacejka.p.D.y(3).*Inc.^2).*Fz;

                Kya = obj.Pacejka.p.K.y(1) .* obj.Pacejka.Fzo .* ( 1 + obj.Pacejka.p.P.y(1).*dPi ) .* ...
                    ( 1 - obj.Pacejka.p.K.y(3).*abs(Inc) ) .* sin( obj.Pacejka.p.K.y(4) .* ...
                    atan( (Fz./obj.Pacejka.Fzo) ./ ...
                    ( ( obj.Pacejka.p.K.y(2) + obj.Pacejka.p.K.y(5).*Inc.^2 ) .* ...
                    ( 1 + obj.Pacejka.p.P.y(2).*dPi ) ) ) );

                Kyg0 = Fz.*(obj.Pacejka.p.K.y(6) + obj.Pacejka.p.K.y(7).*dFz ) .* ...
                    (1 + obj.Pacejka.p.P.y(5).*dPi );

                By = Kya ./ ( Cy.*Dy );

                Vyg = Fz.*(obj.Pacejka.p.V.y(3) + obj.Pacejka.p.V.y(4).*dFz ).*Inc;

                Vy = Fz.*(obj.Pacejka.p.V.y(1) + obj.Pacejka.p.V.y(2).*dFz ) + Vyg;

                Hy = (obj.Pacejka.p.H.y(1) + obj.Pacejka.p.H.y(2).*dFz ) .* (Kyg0.*Inc - Vyg ) ./ Kya;

                Ey = ( obj.Pacejka.p.E.y(1) + obj.Pacejka.p.E.y(2).*dFz ) .* ...
                    ( 1 + obj.Pacejka.p.E.y(5).*Inc.^2 - ...
                    ( obj.Pacejka.p.E.y(3) + obj.Pacejka.p.E.y(4).*Inc ).*sign(Alpha + Hy) );

                % Evaluate Lateral Force
                switch Fidelity.Pure
                    case 'Linear'
                        Fy0 = Kya .* Alpha;
                    case 'Pacejka'
                        Fy0 = (Dy .* sin( Cy .* atan( (1-Ey) .* By.*( Alpha + Hy ) + ...
                            Ey .* atan( By.*( Alpha + Hy ) ) ) ) + Vy) .* (-1).^(mod(i+1,2));
                end

                % Evaluate Null Slip Angle
                if strcmpi( Fidelity.Combined, 'MNC' )
                    Alpha0 = -Vy ./Kya - Hy;
%                  Opts = optimoptions( 'fsolve', 'Display', 'off', 'Algorithm', 'Levenberg-Marquardt' );
%                    Alpha0 = fsolve( @(a) (Dy .* sin( Cy .* atan( (1-Ey) .* By.*(a +Hy) + ...
%                           Ey.*atan( By.*(a + Hy) ) ) ) + Vy) .* (-1).^(mod(i+1,2)), zeros(size(Dy)), Opts );
                      
                else
                    Alpha0 = [];
                end
            end

            %%% Evaluate Pure Slip Aligning Moment (Mz0)
            function [Mz0, Bt, Ct, Dt, Et, Ht, Br, Cr, Dr, Hf] = ...
                    EvaluateMz0( obj, Alpha, ~, Fz, dFz, ~, dPi, Inc, ~, i, ~, ...
                        Fyo, By, Cy, Kya, Hy, Vy )
                % Evaluate P6 Pacejka
                Bt = (obj.Pacejka.q.B.z(1) + obj.Pacejka.q.B.z(2).*dFz + ...
                    obj.Pacejka.q.B.z(3).*dFz.^2) .* ...
                    (1 + obj.Pacejka.q.B.z(5).*abs(Inc) + ...
                    obj.Pacejka.q.B.z(6).*Inc.^2);

                Ct = obj.Pacejka.q.C.z(1);

                Dt = (obj.Pacejka.Ro ./ obj.Pacejka.Fzo) .* ...
                    ( obj.Pacejka.q.D.z(5) .* (obj.Pacejka.Fzo ./ obj.Pacejka.Ro) + ...
                    obj.Pacejka.q.D.z(1).*Fz + obj.Pacejka.q.D.z(2).*Fz.*dFz ) .* ...
                    ( 1 - obj.Pacejka.p.P.z(1).*dPi ) .* ...
                    ( 1 + obj.Pacejka.q.D.z(3).*abs(Inc) + ...
                    obj.Pacejka.q.D.z(4).*Inc.^2 );

                Ht = obj.Pacejka.q.H.z(1) + ...
                    obj.Pacejka.q.H.z(2).*dFz + (obj.Pacejka.q.H.z(3) + ...
                    obj.Pacejka.q.H.z(4).*dFz).*Inc;

                Et = (obj.Pacejka.q.E.z(1) + obj.Pacejka.q.E.z(2).*dFz + ...
                    obj.Pacejka.q.E.z(3).*dFz.^2) .* ...
                    (1 + (obj.Pacejka.q.E.z(4) + obj.Pacejka.q.E.z(5).*Inc).*(2/pi).* ...
                    atan( Bt .* Ct .* (Alpha+Ht) ) );

                Br = obj.Pacejka.q.B.z(10) .* By .* Cy;

                Cr = 1;

                Dr = obj.Pacejka.Ro .* Fz .* ( ( obj.Pacejka.q.D.z(6) + ...
                    obj.Pacejka.q.D.z(7).*dFz ) + ...
                    ( ( obj.Pacejka.q.D.z(8) + obj.Pacejka.q.D.z(9).*dFz ) .* ...
                    ( 1 + obj.Pacejka.p.P.z(2).*dPi ) + ...
                    ( obj.Pacejka.q.D.z(10) + obj.Pacejka.q.D.z(11).*dFz ) .* ...
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
            function Mz = EvaluateMz( obj, Alpha, Kappa, ~, dFz, ~, ~, Inc, ~, i, ~, ...
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
                s = obj.Pacejka.Ro .* (obj.Pacejka.s.s.z(1) + ...
                    obj.Pacejka.s.s.z(2) .* (Fy./obj.Pacejka.Fzo) + ...
                    (obj.Pacejka.s.s.z(3) + obj.Pacejka.s.s.z(4) .* dFz) .* Inc );

                % Evaluate Residual Aligning Moment
                Mzr = Dr .* cos( Cr.*atan( Br.*(Alphareq) ) ) .* cos( Alpha );

                % Evaluate Pure Slip Aligning Moment (Mz0)
                Mz = -t .* Fyp + Mzr + s .* Fx;
            end

            %%% Evaluate Overturning Moment
            function Mx = EvaluateMx( obj, ~, ~, Fz, ~, ~, dPi, Inc, ~, ~, ~, Fy )
%                   Mx = 0;
                
                Mx = obj.Pacejka.Ro .* Fz .* ( obj.Pacejka.q.s.x(1) - ...
                    obj.Pacejka.q.s.x(2) .* Inc .* ( 1 + obj.Pacejka.p.P.Mx(1).*dPi ) + ...
                    obj.Pacejka.q.s.x(3) .* Fy ./ obj.Pacejka.Fzo + ...
                    obj.Pacejka.q.s.x(4) .* cos( obj.Pacejka.q.s.x(5) .* ...
                    atan( (obj.Pacejka.q.s.x(6) .* Fz./obj.Pacejka.Fzo) ).^2 ) .* ...
                    sin( obj.Pacejka.q.s.x(7) .* Inc + obj.Pacejka.q.s.x(8) .* ...
                    atan( obj.Pacejka.q.s.x(9) .* Fy./obj.Pacejka.Fzo ) ) + ...
                    obj.Pacejka.q.s.x(10) .* atan( obj.Pacejka.q.s.x(11) .* Fz./obj.Pacejka.Fzo) .* Inc );
                
            end

            %%% Evaluate Rolling Resistance
            function My = EvaluateMy( obj, Alpha, Kappa, Fz, dFz, Pi, dPi, Inc, V, ~, Fidelity, Fx )
%                 My = 0;
                Re = obj.Radius.Effective( Kappa, Fz, Pi );
                Rl = obj.Radius.Loaded( Fz, Pi );
    
                My = Fx .* (Re - Rl) + Fz .* Re .* ( obj.Pacejka.q.s.y(1) + ...
                    obj.Pacejka.q.s.y(3) .* (V ./ obj.Pacejka.Vo) + ...
                    obj.Pacejka.q.s.y(4) .* (V ./ obj.Pacejka.Vo).^4 );
            end
        end
    end
end