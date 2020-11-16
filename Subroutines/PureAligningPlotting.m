function [] = PureAligningPlotting( Mesh, Raw, Nominal, Tire )
% Plots Results from Pure Aligning (Mzo) MF6.1 Fitting

%% Declare Global Variables
global Figure

%% Evaluate Variant Fit
Mzo = VariantEval( Tire );

%% Operating Case Plots
Figure.Mzo.Cases = figure(); close(gcf);

for p = 1 : size( Nominal, 1 )
    Figure.Mzo.Cases(p) = figure( ...
        'Name', ['Aligning Moment Cases, Pressure (p): ', num2str(Mesh(p,1,1).Pressure), ' [psi]'], ...
        'NumberTitle', 'off');
    
    for z = 1 : size( Nominal, 2 )
        for c = 1 : size( Nominal, 3 )
            subplot( size( Nominal, 3 ), size( Nominal, 2 ), ...
                sub2ind( [size( Nominal, 2 ), size( Nominal, 3 )], z, c ) );
        
            plot( Raw(p,z,c).Slip, Raw(p,z,c).Moment, 'k.' ); hold on;
            
            %fplot( NominalEval( Nominal(p,z,c).C, Nominal(p,z,c).D, ...
            %    Nominal(p,z,c).E, Nominal(p,z,c).K, ...
            %    Nominal(p,z,c).H, Nominal(p,z,c).V ), ...
            %    [-15 15], 'g-.' )
            fplot( @(Slip) Mzo( Mesh(p,z,c).Pressure, Mesh(p,z,c).Load, ...
                Mesh(p,z,c).Camber, Slip), [-15 15], 'g' )
            
            plot( Raw(p,z,c).Slip, Nominal(p,z,c).Residual, 'r.')
            plot( Raw(p,z,c).Slip, ...
                Raw(p,z,c).Moment-Mzo( Mesh(p,z,c).Pressure, Mesh(p,z,c).Load, Mesh(p,z,c).Camber, Raw(p,z,c).Slip), 'y.')
            
            xlabel( 'Slip Angle ($\alpha$) [$deg$]' )
            ylabel( 'Aligning Moment ($M_{z}$) [$N-m$]' )
            title( { ['Normal Load ($F_{z}$): $', num2str(round(Mesh(p,z,c).Load,1)), '$ [$N$]'], ...
                ['Camber ($\gamma$): $', num2str(Mesh(p,z,c).Camber), '$ [$deg$]'] } )
            
            if all([p z c] == ones(3,1))
                legend( {'Raw Data', 'Nominal Fit', 'Variant Fit',...
                    'Nominal Residual', 'Variant Residual'} )
            end
        end
    end
    
    sgtitle( {'Nominal P6 Pacejka Fits', ...
        ['Pressure ($P_{i}$): $', num2str(Mesh(p,1,1).Pressure), '$ [$psi$]'] } )
    
    Figure.Mzo.Cases(p).WindowState = 'minimized';
end

%% Surface Plotting
Figure.Mzo.Surfaces = figure( 'Name', 'Aligning Moment Surfaces', 'NumberTitle', 'off');
    
for p = 1 : size( Raw, 1 )
    for c = 1 : size( Raw, 3 )
        subplot( size( Raw, 3 ), size( Raw, 1 ), ...
            sub2ind( [size( Raw, 1 ), size( Raw, 3 )], p, c ) );
        
        plot3( [Raw(p,:,c).Load], [Raw(p,:,c).Slip], [Raw(p,:,c).Moment], 'k.' ); hold on;
        fsurf( @(Fz, Slip) Mzo( Mesh(p,1,c).Pressure, Fz, ...
            Mesh(p,1,c).Camber, Slip ), [0 2500 -15 15] )
        
        xlabel( 'Normal Load ($F_{z}$) [$N$]' )
        ylabel( 'Slip Angle ($\alpha$) [$deg$]' )
        zlabel( 'Aligning Moment ($M_{z}$) [$N-m$]' )
        title( { ['Pressure ($P_{i}$): $', num2str(Mesh(p,1,c).Pressure), '$ [$psi$]'], ...
                ['Camber ($\gamma$): $', num2str(Mesh(p,1,c).Camber), '$ [$deg$]'] } )
    end
end

sgtitle( 'Pure Aligning MF6.1 Pacejka Fit' )
Figure.Mzo.Surfaces.WindowState = 'minimized';

%% Local Functions
    function Mzo = NominalEval( C, D, E, K, H, V )
        % Magic Formula Evaluation
        B = K ./ ( C.*D + 0.1);
        
        Mzo = @(Slip) ( D.*sin( C.*atan( ...
            ( 1 - E ).*( B.*( Slip + H ) ) + ...
            E.*atan( B .* ( Slip + H ) ) ) ) ) + V;
    end

    function Mzo = VariantEval( Tire )
        % Operating Condition Functions
        dFz = @(Fz) (Fz - Tire.Fzo) ./ Tire.Fzo;
        dPi = @(Pi) (Pi - Tire.Pio) ./ Tire.Pio;
        
        % Lateral Force Evaluation
        Cy = Tire.p.C.y(1);
        
        Dy = @(Pi, Fz, Gam) (Tire.p.D.y(1) + Tire.p.D.y(2).*dFz(Fz)) .* ...
            (1 + Tire.p.P.y(3).*dPi(Pi) + Tire.p.P.y(4).*dPi(Pi).^2) .* ...
            (1 - Tire.p.D.y(3).*Gam.^2).*Fz;
        
        Kya = @(Pi, Fz, Gam) Tire.p.K.y(1) .* Tire.Fzo .* ( 1 + Tire.p.P.y(1).*dPi(Pi) ) .* ...
            ( 1 - Tire.p.K.y(3).*abs(Gam) ) .* sin( Tire.p.K.y(4) .* ...
            atan( (Fz./Tire.Fzo) ./ ...
            ( ( Tire.p.K.y(2) + Tire.p.K.y(5).*Gam.^2 ) .* ( 1 + Tire.p.P.y(2).*dPi(Pi) ) ) ) );
        
        Kyg0 = @(Pi, Fz) Fz.*(Tire.p.K.y(6) + Tire.p.K.y(7).*dFz(Fz) ) .* (1 + Tire.p.P.y(5).*dPi(Pi) );
        
        By = @(Pi, Fz, Gam) Kya(Pi, Fz, 0) ./ ( Cy.*Dy(Pi, Fz, 0) );
        
        Vyg = @(Fz, Gam) Fz.*(Tire.p.V.y(3) + Tire.p.V.y(4).*dFz(Fz) ).*Gam;
        
        Vy = @(Fz, Gam) Fz.*(Tire.p.V.y(1) + Tire.p.V.y(2).*dFz(Fz) ) + Vyg(Fz, Gam);
        
        Hy = @(Pi, Fz, Gam) (Tire.p.H.y(1) + Tire.p.H.y(2).*dFz(Fz) ) .* ...
            (Kyg0(Pi, Fz).*Gam - Vyg(Fz, Gam) ) ./ Kya(Pi, Fz, Gam);
        
        Ey = @(Fz, Gam, Slip, Hy) ( Tire.p.E.y(1) + Tire.p.E.y(2).*dFz(Fz) ) .* ...
            ( 1 + Tire.p.E.y(5).*Gam.^2 - ...
            ( Tire.p.E.y(3) + Tire.p.E.y(4).*Gam ).*sign(Slip + Hy) );

        Fyo = @(Pi, Fz, Gam, Slip) Dy(Pi, Fz, Gam) .* ...
            sin( Cy .* atan( (1-Ey(Fz, Gam, Slip, Hy(Pi, Fz, Gam) )) .* ...
            By(Pi, Fz, Gam).*(Slip + Hy(Pi, Fz, Gam) ) + ...
            Ey(Fz, Gam, Slip, Hy(Pi, Fz, Gam) ).*atan( ...
            By(Pi, Fz, Gam).*(Slip + Hy(Pi, Fz, Gam) ) ) ) ) + Vy(Fz, Gam);
        
        % Pneumatic Trail & Residual Moment
        Bt = @(Fz, Gam) (Tire.q.B.z(1) + Tire.q.B.z(2).*dFz(Fz) + Tire.q.B.z(3).*dFz(Fz).^2) .* ...
            (1 + Tire.q.B.z(5).*abs(Gam) + Tire.q.B.z(6).*Gam.^2);
        
        Ct = Tire.q.C.z(1);
        
        Dt = @(Pi, Fz, Gam) Tire.Ro .* (Fz./Tire.Fzo) .* ...
            ( Tire.q.D.z(1) + Tire.q.D.z(2).*dFz(Fz) ) .* ( 1 - Tire.p.P.z(1) .*dPi(Pi) ) .* ...
            ( 1 + Tire.q.D.z(3).*abs(Gam) + Tire.q.D.z(4).*Gam.^2 );
        
        Ht = @(Fz, Gam) Tire.q.H.z(1) + Tire.q.H.z(2).*dFz(Fz) + ...
            (Tire.q.H.z(3) + Tire.q.H.z(4).*dFz(Fz)).*Gam;
        
        Et = @(Fz, Gam, Slip) (Tire.q.E.z(1) + Tire.q.E.z(2).*dFz(Fz) + Tire.q.E.z(3).*dFz(Fz).^2) .* ...
            (1 + (Tire.q.E.z(4) + Tire.q.E.z(5).*Gam).*(2/pi).* ...
            atan( Bt(Fz, Gam) .* Ct .* deg2rad(Slip + Ht(Fz, Gam) ) ) );
        
        Br = @(Pi, Fz, Gam) Tire.q.B.z(10) .* By(Pi, Fz, Gam) .* Cy;
        
        Cr = 1;
        
        Dr = @(Pi, Fz, Gam) Tire.Ro .* Fz .* ( ( Tire.q.D.z(6) + Tire.q.D.z(7).*dFz(Fz) ) + ...
            ( ( Tire.q.D.z(8) + Tire.q.D.z(9).*dFz(Fz) ) .* ( 1 + Tire.p.P.z(2).*dPi(Pi) ) + ...
            ( Tire.q.D.z(10) + Tire.q.D.z(11).*dFz(Fz) ) .* abs(Gam) ).*Gam );
        
        Hf = @(Pi, Fz, Gam) Hy(Pi, Fz, Gam) + Vy(Fz, Gam) ./ Kya(Pi, Fz, Gam);
        
        % Evaluate Functions & Error
        t0 = @(Pi, Fz, Gam, Slip) Dt(Pi, Fz, Gam) .* cos( Ct .* ...
            atan( (1-Et(Fz, Gam, Slip)).*(Bt(Fz, Gam).*(Slip+Ht(Fz, Gam)) + ...
            Et(Fz, Gam, Slip).*atan( Bt(Fz, Gam).*(Slip+Ht(Fz,Gam)) ) ) ) ) .* cosd( Slip );
        
        Mzro = @(Pi, Fz, Gam, Slip) Dr(Pi, Fz, Gam) .* cos( Cr .* ...
            atan( Br(Pi, Fz, Gam).*(Slip+Hf(Pi, Fz, Gam)) ) ) .* cosd( Slip );
        
        Mzo = @(Pi, Fz, Gam, Slip) -t0(Pi, Fz, Gam, Slip) .* ...
            Fyo(Pi, Fz, Gam, Slip) + Mzro(Pi, Fz, Gam, Slip);
    end
end