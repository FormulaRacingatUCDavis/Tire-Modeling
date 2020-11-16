classdef FRUCDTire
    properties
        Name
        Date
        Model
        Source
        Parameters
    end
​
    methods
        function obj = FRUCDTire(Name, Model, Source, Parameters)
            if nargin > 0
                obj.Name = Name;
                obj.Date = datestr(now,'yyyy-mm-dd_HH_MM');
                obj.Model = Model;
                obj.Source = Source;
                obj.Parameters = Parameters;
            end
        end
        
        function [Fx, Fy, Mz, Mx] = ...
                EvalMF6_1( obj, Alpha, Kappa, Fz, Pi, Gam, Side, Comb )
            
            Alpha = deg2rad(Alpha);
            
            [Fx0, Kxk, Kappa0] = Fx0MF6_1( obj, Alpha, Kappa, Fz, Pi, Gam );
            
            [Fy0, By, Cy, Kya, Hy, Vy, Alpha0] = Fy0MF6_1( obj, Alpha, Kappa, Fz, Pi, Gam, Side );
            
            switch Comb
                case 'Pure'
                    Fx = Fx0;
                    Fy = Fy0;
                    
                case 'MNC'
                    Fx = abs(Fx0 .* Fy0 ./ sqrt( (Kappa - Kappa0).^2 .* Fy0.^2 + Fx0.^2 .* tan(Alpha - Alpha0).^2 ) .* ...
                        sqrt( (Kappa - Kappa0).^2 .* Kya.^2 + (1 - abs(Kappa - Kappa0)).^2 .* cos(Alpha - Alpha0).^2 .* Fx0.^2 ) ./ ...
                        Kya) .* sign(Fx0);
                    
                    Fy = abs( Fx0 .* Fy0 ./ sqrt( (Kappa - Kappa0).^2 .* Fy0.^2 + Fx0.^2 .* tan(Alpha - Alpha0).^2 ) .* ...
                        sqrt( (1 - abs(Kappa - Kappa0)).^2 .* cos(Alpha - Alpha0).^2 .* Fy0.^2 + sin(Alpha - Alpha0).^2 .* Kxk.^2 ) ./ ...
                        ( Kxk .* cos(Alpha - Alpha0)) ) .* sign(Fy0);
        
                case 'MF6_1'
                    % Not Implemented
            end
            
            [Fyo, ~, ~, ~, ~, ~] = Fy0MF6_1( obj, Alpha, Kappa, Fz, Pi, 0, Side );
            
            Mz = Mz0MF6_1( obj, Alpha, Kappa, Fz, Pi, Gam, Fx, Fyo, ...
                By, Cy, Kya, Hy, Vy );
            
            Mx = MxMF6_1( obj, Alpha, Kappa, Fz, Pi, Gam, Fx, Fy );
        end
        
        function [Fx, Kxk, Kappa0] = Fx0MF6_1( obj, ~, Kappa, Fz, Pi, Gam )
            dFz = (Fz - obj.Parameters.Fzo) ./ obj.Parameters.Fzo;
            dPi = (Pi - obj.Parameters.Pio) ./ obj.Parameters.Pio;
​
            Cx = obj.Parameters.p.C.x(1) .* obj.Parameters.L.C.x;
​
            Dx = (obj.Parameters.p.D.x(1) + obj.Parameters.p.D.x(2).*dFz) .* ...
                (1 + obj.Parameters.p.P.x(3).*dPi + obj.Parameters.p.P.x(4).*dPi.^2) .* ...
                (1 - obj.Parameters.p.D.x(3).*Gam.^2).*Fz .* obj.Parameters.L.mu.x;
​
            Ex = ( obj.Parameters.p.E.x(1) + obj.Parameters.p.E.x(2).*dFz ...
                + obj.Parameters.p.E.x(3).*dFz.^2 ) .* ( 1 - obj.Parameters.p.E.x(4).*sign(Kappa) ) .* ...
                obj.Parameters.L.E.x;
​
            Kxk = Fz.*(obj.Parameters.p.K.x(1) + obj.Parameters.p.K.x(2).*dFz ) .* ...
                exp( obj.Parameters.p.K.x(3) .* dFz ) .* ...
                (1 + obj.Parameters.p.P.x(1).*dPi + obj.Parameters.p.P.x(2).*dPi.^2) .* ...
                obj.Parameters.L.K.x.k;
​
            Bx = Kxk ./ ( Cx.*Dx );
​
            Vx = Fz.*(obj.Parameters.p.V.x(1) + obj.Parameters.p.V.x(2).*dFz ) .* ...
                obj.Parameters.L.V.x;
​
            Hx = (obj.Parameters.p.H.x(1) + obj.Parameters.p.H.x(2).*dFz) .* ...
                obj.Parameters.L.H.x;
​
            Fx = Dx .* sin( Cx .* atan( (1-Ex) .* Bx.*(Kappa + Hx) + ...
                 Ex.*atan( Bx.*(Kappa + Hx) ) ) ) + Vx;
           
            Opts = optimoptions( 'fsolve', 'Display', 'off' );
            Kappa0 = fsolve( @(k) Dx .* sin( Cx .* atan( (1-Ex) .* Bx.*(k + Hx) + ...
                 Ex.*atan( Bx.*(k + Hx) ) ) ) + Vx, zeros(size(Dx)), Opts );
        end
        
        function [Fy, By, Cy, Kya, Hy, Vy, Alpha0] = Fy0MF6_1( obj, Alpha, ~, Fz, Pi, Gam, Side )
            dFz = (Fz - obj.Parameters.Fzo) ./ obj.Parameters.Fzo;
            dPi = (Pi - obj.Parameters.Pio) ./ obj.Parameters.Pio;
​
            Cy = obj.Parameters.p.C.y(1) .* obj.Parameters.L.C.y;
​
            Dy = (obj.Parameters.p.D.y(1) + obj.Parameters.p.D.y(2).*dFz) .* ...
                (1 + obj.Parameters.p.P.y(3).*dPi + obj.Parameters.p.P.y(4).*dPi.^2) .* ...
                (1 - obj.Parameters.p.D.y(3).*Gam.^2).*Fz .* obj.Parameters.L.mu.y;
​
            Kya = obj.Parameters.p.K.y(1) .* obj.Parameters.Fzo .* ( 1 + obj.Parameters.p.P.y(1).*dPi ) .* ...
                ( 1 - obj.Parameters.p.K.y(3).*abs(Gam) ) .* sin( obj.Parameters.p.K.y(4) .* ...
                atan( (Fz./obj.Parameters.Fzo) ./ ...
                ( ( obj.Parameters.p.K.y(2) + obj.Parameters.p.K.y(5).*Gam.^2 ) .* ...
                ( 1 + obj.Parameters.p.P.y(2).*dPi ) ) ) ) .* obj.Parameters.L.K.y.a;
​
            Kyg0 = Fz.*(obj.Parameters.p.K.y(6) + obj.Parameters.p.K.y(7).*dFz ) .* ...
                (1 + obj.Parameters.p.P.y(5).*dPi ) .* obj.Parameters.L.K.y.g;
​
            By = Kya ./ ( Cy.*Dy );
​
            Vyg = Fz.*(obj.Parameters.p.V.y(3) + obj.Parameters.p.V.y(4).*dFz ).*Gam .* ...
                obj.Parameters.L.K.y.g .* obj.Parameters.L.mu.y;
​
            Vy = Fz.*(obj.Parameters.p.V.y(1) + obj.Parameters.p.V.y(2).*dFz ) .* ...
                obj.Parameters.L.V.y .* obj.Parameters.L.mu.y + Vyg;
​
            Hy = (obj.Parameters.p.H.y(1) + obj.Parameters.p.H.y(2).*dFz ) .* ...
                obj.Parameters.L.H.y + (Kyg0.*Gam - Vyg ) ./ Kya;
​
            Ey = ( obj.Parameters.p.E.y(1) + obj.Parameters.p.E.y(2).*dFz ) .* ...
                ( 1 + obj.Parameters.p.E.y(5).*Gam.^2 - ...
                ( obj.Parameters.p.E.y(3) + obj.Parameters.p.E.y(4).*Gam ).*sign(Alpha + Hy) ) .* ...
                obj.Parameters.L.E.y;
            
            if Side == 'L'
                Fy = Dy .* sin( Cy .* atan( (1-Ey) .* By.*(Alpha + Hy ) + ...
                    Ey.*atan( By.*(Alpha + Hy ) ) ) ) + Vy;
​
                Opts = optimoptions( 'fsolve', 'Display', 'off' );
                Alpha0 = fsolve( @(a) Dy .* sin( Cy .* atan( (1-Ey) .* By.*(a + Hy ) + ...
                    Ey.*atan( By.*(a + Hy ) ) ) ) + Vy, zeros(size(Dy)), Opts );
            else
                Fy = Dy .* sin( Cy .* atan( (1-Ey) .* By.*(Alpha - Hy ) + ...
                    Ey.*atan( By.*(Alpha - Hy ) ) ) ) - Vy;
            
                Opts = optimoptions( 'fsolve', 'Display', 'off' );
                Alpha0 = fsolve( @(a) Dy .* sin( Cy .* atan( (1-Ey) .* By.*(a - Hy ) + ...
                    Ey.*atan( By.*(a - Hy ) ) ) ) - Vy, zeros(size(Dy)), Opts );
            end
        end
        
        function Mz = Mz0MF6_1( obj, Alpha, ~, Fz, Pi, Gam, ~, Fy, ...
                By, Cy, Kya, Hy, Vy )
            dFz = (Fz - obj.Parameters.Fzo) ./ obj.Parameters.Fzo;
            dPi = (Pi - obj.Parameters.Pio) ./ obj.Parameters.Pio;
            
            Bt = (obj.Parameters.q.B.z(1) + obj.Parameters.q.B.z(2).*dFz + ...
                obj.Parameters.q.B.z(3).*dFz.^2) .* ...
                (1 + obj.Parameters.q.B.z(5).*abs(Gam) + ...
                obj.Parameters.q.B.z(6).*Gam.^2);
        
            Ct = obj.Parameters.q.C.z(1);
​
            Dt = (obj.Parameters.Ro ./ obj.Parameters.Fzo) .* ...
                ( obj.Parameters.q.D.z(5) .* (obj.Parameters.Fzo ./ obj.Parameters.Ro) + ...
                obj.Parameters.q.D.z(1).*Fz + obj.Parameters.q.D.z(2).*Fz.*dFz ) .* ...
                ( 1 - obj.Parameters.p.P.z(1).*dPi ) .* ...
                ( 1 + obj.Parameters.q.D.z(3).*abs(Gam) + ...
                obj.Parameters.q.D.z(4).*Gam.^2 );
​
            Ht = obj.Parameters.q.H.z(1) + ...
                obj.Parameters.q.H.z(2).*dFz + (obj.Parameters.q.H.z(3) + ...
                obj.Parameters.q.H.z(4).*dFz).*Gam;
​
            Et = (obj.Parameters.q.E.z(1) + obj.Parameters.q.E.z(2).*dFz + ...
                obj.Parameters.q.E.z(3).*dFz.^2) .* ...
                (1 + (obj.Parameters.q.E.z(4) + obj.Parameters.q.E.z(5).*Gam).*(2/pi).* ...
                atan( Bt .* Ct .* (Alpha+Ht) ) );
​
            Br = obj.Parameters.q.B.z(10) .* By .* Cy;
​
            Cr = 1;
​
            Dr = obj.Parameters.Ro .* Fz .* ( ( obj.Parameters.q.D.z(6) + ...
                obj.Parameters.q.D.z(7).*dFz ) + ...
                ( ( obj.Parameters.q.D.z(8) + obj.Parameters.q.D.z(9).*dFz ) .* ...
                ( 1 + obj.Parameters.p.P.z(2).*dPi ) + ...
                ( obj.Parameters.q.D.z(10) + obj.Parameters.q.D.z(11).*dFz ) .* ...
                abs(Gam) ).*Gam );
​
            Hf = Hy + Vy ./ Kya;
​
            % Evaluate Functions & Error
            t0 = Dt.*cos( Ct.*atan( (1-Et).*(Bt.*(Alpha+Ht) + ...
                Et.*atan( Bt.*(Alpha+Ht) ) ) ) ) .* cos( Alpha );
​
            Mzro = Dr .* cos( Cr.*atan( Br.*(Alpha+Hf) ) ) .* cos( Alpha );
​
            Mz = -t0 .* Fy + Mzro;
        end
        
        function Mx = MxMF6_1( ~, Alpha, Kappa, Fz, Pi, Gam, Fx, Fy )
           Mx = 0 * Alpha .* Kappa .* Fz .* Pi .* Gam .* Fy .* Fx;
           
           %{
           Mxo = @(Pi, Fz, Gam, Slip) obj.Parameters.Ro .* Fz .* ( obj.Parameters.q.s.x(1) - ...
           Tire.q.s.x(2) .* Gam .* ( 1 + Tire.p.P.Mx(1).*dPi(Pi) ) + ...
           Tire.q.s.x(3) .* Fyo(Pi, Fz, Gam, Slip)./Tire.Fzo + ...
           Tire.q.s.x(4) .* cos( Tire.q.s.x(5) .* atan( (Tire.q.s.x(6) .* Fz./Tire.Fzo).^2 ) ) .* ...
           sin( Tire.q.s.x(7) .* Gam + Tire.q.s.x(8) .* atan( Tire.q.s.x(9) .* Fyo( Pi, Fz, Gam, Slip)./Tire.Fzo ) ) + ...
           Tire.q.s.x(10) .* atan( Tire.q.s.x(11) .* Fz./Tire.Fzo) .* Gam );
           %}
        end
    end
end