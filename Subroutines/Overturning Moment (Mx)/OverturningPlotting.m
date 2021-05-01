function [] = OverturningPlotting(Mesh, Raw, Variant, Tire, Figure )
%% Overturning Plotting = Plots Results from Overturning (Mx) Fitting
% Plots Variant Values from fitting process for Overturning Process for the
% equations given in Pacejka's "Tire and Vehicle Dynamics" [3rd Edition] in
% section 4.3.2 (page 176). 

% Inputs 
%   Mesh    - Operating Condition Space
%   Raw     - Allocated Data
%   Variant - Parameter Fit Values for Pacejka Coefficients
%   Tire    - Tire Model
%   Figure  - Stores Model Figures

% Outputs
%  Surface Plots of Overturning Fitting

% Authors:
% Carlos Lopez (calopez@ucdavis.edu) [Dec 2020 - June 2022]

% Last Updated: 17- APR - 2021
%% Evaluate Variant Surface
[Mx] = VariantEval( Tire );



%% Variant Surface Plotting
Figure.Mx.Surfaces = figure( 'Name'       , 'Overturning Moment Surfaces', ...
                              'NumberTitle', 'off', ...
                              'Visible'    , 'on' );

for p = 1 : size( Raw, 1 )
    for c = 1 : size( Raw, 3 )
        subplot( size( Raw, 3 ), size( Raw, 1 ), ...
            sub2ind( [size( Raw, 1 ), size( Raw, 3 )], p, c ) );
        
        plot3( [Raw(p,:,c).Load], rad2deg([Raw(p,:,c).Alpha]), [Raw(p,:,c).Moment], 'k.' ); hold on;
        fsurf( @(Fz, Alpha) Mx(Alpha,0, Fz, Mesh(p,1,c).Pressure, ...
            Mesh(p,1,c).Inclination), [0 2500 -15 15] )
    
        xlabel( 'Normal Load ($F_{z}$) [$N$]' )
        ylabel( 'Slip Angle ($\alpha$) [$deg$]' )
        zlabel( 'Overturning Moment ($M_{x}$) [$Nm$]' )
        title( { ['Pressure ($P_{i}$): $'    , num2str(Mesh(p,1,c).Pressure)   , '$ [$kPa$]'], ...
                 ['Inclination ($\gamma$): $', num2str(Mesh(p,1,c).Inclination), '$ [$deg$]'] } )
 
    end
end

sgtitle( 'Overturning MF6.1 Pacejka Fit' )
Figure.Mx.Surfaces.WindowState = Figure.State;
%% Local Functions

function [Mx] = VariantEval( Tire )
    
        dPi = @(Pi) (Pi - Tire.Pacejka.Pio) ./ Tire.Pacejka.Pio;
        dFz = @(Fz) (Fz - Tire.Pacejka.Fzo) ./ Tire.Pacejka.Fzo;
        
        % Lateral Force Evaluation (Inclination = 0)
        Cy = Tire.Pacejka.p.C.y(1);
        
        Dy = @(Fz, Pi, Gam) (Tire.Pacejka.p.D.y(1) + Tire.Pacejka.p.D.y(2).*dFz(Fz)) .* ...
            (1 + Tire.Pacejka.p.P.y(3).*dPi(Pi) + Tire.Pacejka.p.P.y(4).*dPi(Pi).^2) .* ...
            (1 - Tire.Pacejka.p.D.y(3).*Gam.^2).*Fz;
        
        Kya = @(Fz, Pi, Gam) Tire.Pacejka.p.K.y(1) .* Tire.Pacejka.Fzo .* ( 1 + Tire.Pacejka.p.P.y(1).*dPi(Pi) ) .* ...
            ( 1 - Tire.Pacejka.p.K.y(3).*abs(Gam) ) .* sin( Tire.Pacejka.p.K.y(4) .* ...
            atan( (Fz./Tire.Pacejka.Fzo) ./ ...
            ( ( Tire.Pacejka.p.K.y(2) + Tire.Pacejka.p.K.y(5).*Gam.^2 ) .* ( 1 + Tire.Pacejka.p.P.y(2).*dPi(Pi) ) ) ) );
        
        Kyg0 = @(Fz,Pi) Fz.*(Tire.Pacejka.p.K.y(6) + Tire.Pacejka.p.K.y(7).*dFz(Fz) ) .* (1 + Tire.Pacejka.p.P.y(5).*dPi(Pi) );
        
        By = @(Fz, Pi, Gam) Kya(Pi, Fz, Gam) ./ ( Cy.*Dy(Fz, Pi, Gam) );
        
        Vyg = @(Fz, Gam) Fz.*(Tire.Pacejka.p.V.y(3) + Tire.Pacejka.p.V.y(4).*dFz(Fz) ).*Gam;
        
        Vy = @(Fz, Gam) Fz.*(Tire.Pacejka.p.V.y(1) + Tire.Pacejka.p.V.y(2).*dFz(Fz) ) + Vyg(Fz, Gam);
        
        Hy = @(Fz, Pi, Gam) (Tire.Pacejka.p.H.y(1) + Tire.Pacejka.p.H.y(2).*dFz(Fz) ) .* ...
            (Kyg0(Fz, Pi).*Gam - Vyg(Fz, Gam) ) ./ Kya(Fz, Pi, Gam);
        
        Ey = @(Alpha, Fz, Gam, Hy) ( Tire.Pacejka.p.E.y(1) + Tire.Pacejka.p.E.y(2).*dFz(Fz) ) .* ...
            ( 1 + Tire.Pacejka.p.E.y(5).*Gam.^2 - ...
            ( Tire.Pacejka.p.E.y(3) + Tire.Pacejka.p.E.y(4).*Gam ).*sign(Alpha + Hy) );

        Fy0 = @(Alpha, Fz, Pi, Gam) Dy(Fz, Pi, Gam) .* ...
            sin( Cy .* atan( (1-Ey(Alpha, Fz, Gam, Hy(Fz, Pi, Gam) )) .* ...
            By(Fz, Pi, Gam).*(Alpha + Hy(Fz, Pi, Gam) ) + ...
            Ey(Alpha, Fz, Gam, Hy(Fz, Pi, Gam) ).*atan( ...
            By(Fz, Pi, Gam).*(Alpha + Hy(Fz, Pi, Gam) ) ) ) ) + Vy(Fz, Gam); 
        
        % Evaluate P6 Pacejka
        Cx = Tire.Pacejka.p.C.x(1) .* Tire.Pacejka.L.C.x;
        
        Dx = @(Fz, Pi, Gam)(Tire.Pacejka.p.D.x(1) + Tire.Pacejka.p.D.x(2).*dFz(Fz)) .* ...
            (1 + Tire.Pacejka.p.P.x(3).*dPi(Pi) + Tire.Pacejka.p.P.x(4).*dPi(Pi).^2) .* ...
            (1 - Tire.Pacejka.p.D.x(3).*Gam.^2).*Fz .* Tire.Pacejka.L.mu.x;
        
        Ex = @(Kappa, Fz) ( Tire.Pacejka.p.E.x(1) + Tire.Pacejka.p.E.x(2).*dFz(Fz) ...
            + Tire.Pacejka.p.E.x(3).*dFz(Fz).^2 ) .* ( 1 - Tire.Pacejka.p.E.x(4).*sign(Kappa) ) .* ...
            Tire.Pacejka.L.E.x;
        
        Kxk = @(Fz, Pi) Fz.*(Tire.Pacejka.p.K.x(1) + Tire.Pacejka.p.K.x(2).*dFz(Fz) ) .* ...
            exp( Tire.Pacejka.p.K.x(3) .* dFz(Fz) ) .* ...
            (1 + Tire.Pacejka.p.P.x(1).*dPi(Pi) + Tire.Pacejka.p.P.x(2).*dPi(Pi).^2) .* ...
            Tire.Pacejka.L.K.x.k;
        
        Bx = @(Fz, Pi, Gam) Kxk(Fz,Pi) ./ ( Cx.*Dx(Fz,Pi, Gam ));
        
        Vx = @(Fz) Fz.*(Tire.Pacejka.p.V.x(1) + Tire.Pacejka.p.V.x(2).*dFz(Fz) ) .* ...
            Tire.Pacejka.L.V.x;
        
        Hx = @(Fz)(Tire.Pacejka.p.H.x(1) + Tire.Pacejka.p.H.x(2).*dFz(Fz)) .* ...
            Tire.Pacejka.L.H.x;        
                
        Fx0 = @(Kappa, Fz, Pi, Gam) Dx(Fz,Pi, Gam) .* sin( Cx .* ... 
            atan( (1-Ex(Kappa, Fz)) .* Bx(Fz,Pi,Gam).*(Kappa + Hx(Fz)) + ...
            Ex(Kappa, Fz).*atan( Bx(Fz, Pi, Gam).*(Kappa + Hx(Fz)) ) ) ) + Vx(Fz); 

        %Defining Alpha0 and Kappa0
        Kappa0 = @(Fz, Pi) -Vx(Fz) ./ Kxk(Fz,Pi) - Hx(Fz);
   
        Alpha0 = @(Fz, Pi, Gam) -Vy(Fz,Gam) ./Kya(Pi, Fz, Gam) - Hy(Pi, Fz, Gam);
        
        % Defining Fy
        Fy = @(Alpha, Kappa, Fz, Pi, Gam ) abs( Fx0(Kappa, Fz, Pi, Gam) .* Fy0(Alpha, Fz, Pi, Gam) ./ ...
            sqrt( (Kappa - Kappa0(Fz, Pi)).^2 .* Fy0(Alpha, Fz, Pi, Gam).^2 + ...
            Fx0(Kappa, Fz, Pi, Gam).^2 .* tan(Alpha - Alpha0(Fz,Pi, Gam)).^2 ) .* ...
            sqrt( (1 - abs(Kappa - Kappa0(Fz,Pi))).^2 .* cos(Alpha - Alpha0(Fz, Pi, Gam)).^2 .*...
            Fy0(Alpha, Fz, Pi, Gam).^2 + sin(Alpha - Alpha0(Fz, Pi, Gam)).^2 .* Kxk(Fz,Pi).^2 ) ./ ...
            ( Kxk(Fz,Pi) .* cos(Alpha - Alpha0(Fz,Pi, Gam))) ) .* sign(Fy0(Alpha, Fz, Pi, Gam));           
        
        Mx = @( Alpha, Kappa,Fz, Pi, Gam) (Tire.Pacejka.Ro .* Tire.Pacejka.Fzo) .* ... 
            ( Tire.Pacejka.q.s.x(1) - ((Tire.Pacejka.q.s.x(2) .* Gam) .* ...
            (1 + Tire.Pacejka.p.P.Mx .* dPi(Pi))) +(Tire.Pacejka.q.s.x(3) .* ...
            (Fy(Alpha, Kappa, Fz, Pi, Gam)./Tire.Pacejka.Fzo)) + ...
            (Tire.Pacejka.q.s.x(4) .* cos(Tire.Pacejka.q.s.x(5) .* ...
            atan(Tire.Pacejka.q.s.x(6) .* (Fz./Tire.Pacejka.Fzo)).^2).* ...
            sin((Tire.Pacejka.q.s.x(7) .* Gam) + ...
            (Tire.Pacejka.q.s.x(8) .* atan( Tire.Pacejka.q.s.x(9) .* ...
            (Fy(Alpha, Kappa, Fz, Pi, Gam)./Tire.Pacejka.Fzo))))) + ... 
            (Tire.Pacejka.q.s.x(10) .* atan(Tire.Pacejka.q.s.x(11) .* ...
            (Fz./Tire.Pacejka.Fzo)) .* Gam));
end
end