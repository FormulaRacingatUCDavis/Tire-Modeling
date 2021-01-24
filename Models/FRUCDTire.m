%% FRUCD Wheel & Tire Class Definition
% Defines wheel & tire parameters for use in vehicle simulation and
% provides methods for evaluating certain properties.
%
% Main Methods:
% - FRUCDTire: Instance Generation
% - ContactPatchLoads: Evaluates Model for Fx, Fy, Mz, Mx, and My
% - RadialDeflection: Evaluates Model for Re, Rl
classdef FRUCDTire
    properties
        Name
        Date
        Source
        Notes
        Pacejka
        Radius
        Mass
        Thermal
    end
    methods
        function obj = FRUCDTire(Name, Source, Notes, Pacejka, Radius, Mass, Thermal)
            if nargin > 0
                obj.Name    = Name;                            % Tire Name
                obj.Date    = datestr(now,'yyyy-mm-dd_HH_MM'); % Generation Date & Time
                obj.Source  = Source;                          % Data Sources     
                obj.Notes   = Notes;                           % Notes About Model
                
                obj.Pacejka = Pacejka;                         % Pacejka Parameters
                obj.Radius  = Radius;                          % Radius Parameters
                obj.Mass    = Mass;                            % Mass Parameters
                obj.Thermal = Thermal;                         % Thermal Parameters
            end
        end    
        
        %%% Evaluate Tire Model for Fx, Fy, Mz, Mx, My
        function [Fx, Fy, Mz, Mx, My] = ContactPatchLoads( obj, Alpha, Kappa, Fz, Pi, Inc, V, i, Fidelity )   
            %%% Slip Angle Conversion
            Alpha = deg2rad(Alpha);
            
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
        end
        
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
                    Fx0 = Kxk * Kappa;
                case 'Pacejka'
                    Fx0 = Dx .* sin( Cx .* atan( (1-Ex) .* Bx.*(Kappa + Hx) + ...
                        Ex.*atan( Bx.*(Kappa + Hx) ) ) ) + Vx;
            end
            
            % Evaluate Null Slip Ratio
            if strcmpi( Fidelity.Combined, 'MNC' )
                Opts = optimoptions( 'fsolve', 'Display', 'off', 'Algorithm', 'Levenberg-Marquardt' );
                Kappa0 = fsolve( @(k) Dx .* sin( Cx .* atan( (1-Ex) .* Bx.*(k + Hx) + ...
                     Ex.*atan( Bx.*(k + Hx) ) ) ) + Vx, zeros(size(Dx)), Opts );
            else
                Kappa0 = [];
            end
        end
        
        %%% Pure Slip Lateral Force Evaluation (Fy0)
        function [Fy0, By, Cy, Kya, Hy, Vy, Alpha0] = ...
                EvaluateFy0( obj, Alpha, ~, Fz, dFz, ~, dPi, Inc, ~, i, Fidelity )            
            % Evaluate P6 Pacejka
            Cy = obj.Pacejka.p.C.y(1) .* obj.Pacejka.L.C.y;
            
            Dy = (obj.Pacejka.p.D.y(1) + obj.Pacejka.p.D.y(2).*dFz) .* ...
                (1 + obj.Pacejka.p.P.y(3).*dPi + obj.Pacejka.p.P.y(4).*dPi.^2) .* ...
                (1 - obj.Pacejka.p.D.y(3).*Inc.^2).*Fz .* obj.Pacejka.L.mu.y;
            
            Kya = obj.Pacejka.p.K.y(1) .* obj.Pacejka.Fzo .* ( 1 + obj.Pacejka.p.P.y(1).*dPi ) .* ...
                ( 1 - obj.Pacejka.p.K.y(3).*abs(Inc) ) .* sin( obj.Pacejka.p.K.y(4) .* ...
                atan( (Fz./obj.Pacejka.Fzo) ./ ...
                ( ( obj.Pacejka.p.K.y(2) + obj.Pacejka.p.K.y(5).*Inc.^2 ) .* ...
                ( 1 + obj.Pacejka.p.P.y(2).*dPi ) ) ) ) .* obj.Pacejka.L.K.y.a;
            
            Kyg0 = Fz.*(obj.Pacejka.p.K.y(6) + obj.Pacejka.p.K.y(7).*dFz ) .* ...
                (1 + obj.Pacejka.p.P.y(5).*dPi ) .* obj.Pacejka.L.K.y.g;
            
            By = Kya ./ ( Cy.*Dy );
            
            Vyg = Fz.*(obj.Pacejka.p.V.y(3) + obj.Pacejka.p.V.y(4).*dFz ).*Inc .* ...
                obj.Pacejka.L.K.y.g .* obj.Pacejka.L.mu.y;
            
            Vy = Fz.*(obj.Pacejka.p.V.y(1) + obj.Pacejka.p.V.y(2).*dFz ) .* ...
                obj.Pacejka.L.V.y .* obj.Pacejka.L.mu.y + Vyg;
            
            Hy = (obj.Pacejka.p.H.y(1) + obj.Pacejka.p.H.y(2).*dFz ) .* ...
                obj.Pacejka.L.H.y + (Kyg0.*Inc - Vyg ) ./ Kya;
            
            Ey = ( obj.Pacejka.p.E.y(1) + obj.Pacejka.p.E.y(2).*dFz ) .* ...
                ( 1 + obj.Pacejka.p.E.y(5).*Inc.^2 - ...
                ( obj.Pacejka.p.E.y(3) + obj.Pacejka.p.E.y(4).*Inc ).*sign(Alpha + Hy) ) .* ...
                obj.Pacejka.L.E.y;
            
            % Evaluate Lateral Force
            switch Fidelity.Pure
                case 'Linear'
                    Fy0 = Kya * Alpha;
                case 'Pacejka'
                    Fy0 = Dy .* sin( Cy .* atan( (1-Ey) .* By.*(Alpha + (-1).^(mod(i+1,2)).*Hy ) + ...
                        Ey.*atan( By.*(Alpha + (-1).^(mod(i+1,2)).*Hy ) ) ) ) + ...
                        (-1).^(mod(i+1,2)).*Vy;
            end

            % Evaluate Null Slip Angle
            if strcmpi( Fidelity.Combined, 'MNC' )
                Opts = optimoptions( 'fsolve', 'Display', 'off', 'Algorithm', 'Levenberg-Marquardt' );
                Alpha0 = fsolve( @(a) Dy .* sin( Cy .* atan( (1-Ey) .* By.*(a + (-1).^(mod(i+1,2)).*Hy ) + ...
                    Ey.*atan( By.*(a + (-1).^(mod(i+1,2)).*Hy ) ) ) ) + (-1).^(mod(i+1,2)).*Vy, zeros(size(Dy)), Opts );
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
                abs(Inc) ).*Inc );
            
            Hf = Hy + Vy ./ Kya;
            
            % Evaluate Trail Function
            t0 = Dt.*cos( Ct.*atan( (1-Et).*(Bt.*(Alpha + (-1).^(mod(i+1,2)).*Ht ) + ...
                Et.*atan( Bt.*(Alpha + (-1).^(mod(i+1,2)).*Ht ) ) ) ) ) .* cos( Alpha );
            
            % Evaluate Residual Aligning Moment
            Mzro = Dr .* cos( Cr.*atan( Br.*(Alpha + (-1).^(mod(i+1,2)).*Hf ) ) ) .* cos( Alpha );
            
            % Evaluate Pure Slip Aligning Moment (Mz0)
            Mz0 = -t0 .* Fyo + Mzro;
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
            Mx = 0;
            %{
            Mx = obj.Pacejka.Ro .* Fz .* ( obj.Pacejka.q.s.x(1) - ...
                obj.Pacejka.q.s.x(2) .* Inc .* ( 1 + obj.Pacejka.p.P.Mx(1).*dPi ) + ...
                obj.Pacejka.q.s.x(3) .* Fy ./ obj.Pacejka.Fzo + ...
                obj.Pacejka.q.s.x(4) .* cos( obj.Pacejka.q.s.x(5) .* ...
                atan( (obj.Pacejka.q.S.x(6) .* Fz./obj.Pacejka.Fzo) ).^2 ) .* ...
                sin( obj.Pacejka.q.s.x(7) .* Inc + obj.Pacejka.q.s.x(8) .* ...
                atan( obj.Pacejka.q.s.x(9) .* Fy./obj.Pacejka.Fzo ) ) + ...
                obj.Pacejka.q.s.x(10) .* atan( obj.Pacejka.q.s.x(11) .* Fz./obj.Pacejka.Fzo) .* Inc );
            %}
        end
    
        %%% Evaluate Rolling Resistance
        function My = EvaluateMy( obj, Alpha, Kappa, Fz, dFz, Pi, dPi, Inc, V, ~, Fidelity, Fx )
            My = 0;
        end
    end
end