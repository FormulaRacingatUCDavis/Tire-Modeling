function [] = PureAligningPlotting( Mesh, Raw, Nominal, ~, ~, Tire )
% Plots Results from Pure Aligning (Mzo) MF6.1 Fitting

%% Declare Global Variables
global Figure

%% Evaluate Variant Fit
[By, Bt, Dt, Et, Ht, Br, Dr, Hf, t0, Mzro, Mzo] = VariantEval( Tire );

%% Operating Case Plots
for p = 1 : size( Nominal, 1 )
    Figure.Mzo.Nominal(p) = figure( ...
        'Name'       , ['Aligning Nominal Fits, Pressure (p): ', num2str(Mesh(p,1,1).Pressure), ' [psi]'], ...
        'NumberTitle', 'off', ...
        'Visible'    , 'on' );
    
    for z = 1 : size( Nominal, 2 )
        for c = 1 : size( Nominal, 3 )
            subplot( size( Nominal, 3 ), size( Nominal, 2 ), ...
                sub2ind( [size( Nominal, 2 ), size( Nominal, 3 )], z, c ) );
        
            plot( rad2deg(Raw(p,z,c).Slip), Raw(p,z,c).Moment, 'k.' ); hold on;
            
            fplot( @(Slip) Mzo( Mesh(p,z,c).Pressure, Mesh(p,z,c).Load, ...
                Mesh(p,z,c).Inclination, Slip), [-15 15], 'g' )
            
            plot( rad2deg(Raw(p,z,c).Slip), Nominal(p,z,c).Residual, 'r.')
            plot( rad2deg(Raw(p,z,c).Slip), Raw(p,z,c).Moment - ...
                Mzo( Mesh(p,z,c).Pressure   , Mesh(p,z,c).Load, ...
                     Mesh(p,z,c).Inclination, Raw(p,z,c).Slip), 'y.')
            
            xlabel( 'Slip Angle ($\alpha$) [$deg$]' )
            ylabel( 'Aligning Moment ($M_{z}$) [$N-m$]' )
            title( { ['Normal Load ($F_{z}$): $', num2str(round(Mesh(p,z,c).Load,1)), '$ [$N$]'], ...
                ['Inclination ($\gamma$): $', num2str(Mesh(p,z,c).Inclination), '$ [$deg$]'] } )
            
            if all([p z c] == ones(3,1))
                legend( {'Raw Data'        , 'Variant Fit'     ,...
                         'Nominal Residual', 'Variant Residual'} )
            end
        end
    end
    
    sgtitle( {'Nominal P6 Pacejka Fits', ...
        ['Pressure ($P_{i}$): $', num2str(Mesh(p,1,1).Pressure), '$ [$psi$]'] } )
    
    Figure.Mzo.Cases(p).WindowState = 'minimized';
end

%% Response Surface Plotting
if strcmpi( Figure.Mode, 'Debug' )
    Figure.Fyo.Response = figure( ...
    'Name'       , 'Aligning Response Surfaces', ...
    'NumberTitle', 'off', ...
    'Visible'    , 'on' );

    ax1 = subplot(20,2,1:2:7);  
    fplot( @(Fz) Bt(Fz, 0), [0 2500], 'm'   ); hold on;
    fplot( @(Fz) Bt(Fz, 2), [0 2500], 'm--' );
    fplot( @(Fz) Bt(Fz, 4), [0 2500], 'm:'  );

    scatter( [Mesh(:).Load], [Nominal(:).Bt], 'k.' ); 
    ylabel( 'Curvature: $B_{t}$' )
    xlabel( 'Normal Load: $F_{z}$' )
         
    ax2 = subplot(20,2,9:2:15);
    fplot( @(Fz) Dt(10, Fz, 0), [0 2500], 'b'   ); hold on;
    fplot( @(Fz) Dt(12, Fz, 0), [0 2500], 'm'   );
    fplot( @(Fz) Dt(12, Fz, 2), [0 2500], 'm--' );
    fplot( @(Fz) Dt(12, Fz, 4), [0 2500], 'm:'  );
    fplot( @(Fz) Dt(14, Fz, 0), [0 2500], 'r'   );

    scatter( [Mesh(:).Load], [Nominal(:).Dt], 'k.' ); 
    ylabel( 'Peak: $D_{t}$' )
    xlabel( 'Normal Load: $F_{z}$' )

    ax3 = subplot(20,2,17:2:23);
    fplot( @(Fz) Et(Fz, 0, 0), [0 2500], 'm'   ); hold on;
    fplot( @(Fz) Et(Fz, 2, 0), [0 2500], 'm--' );
    fplot( @(Fz) Et(Fz, 4, 0), [0 2500], 'm:'  );

    scatter( [Mesh(:).Load], [Nominal(:).Et], 'k.' ); 
    ylabel( 'Asymptote: $E_{t}$' )
    xlabel( 'Normal Load: $F_{z}$' )

    ax4 = subplot(20,2,25:2:31);
    fplot( @(Fz) Ht(Fz, 0), [0 2500], 'm'   ); hold on;
    fplot( @(Fz) Ht(Fz, 2), [0 2500], 'm--' );
    fplot( @(Fz) Ht(Fz, 4), [0 2500], 'm:'  );

    scatter( [Mesh(:).Load], [Nominal(:).Ht], 'k.' ); 
    ylabel( 'Shift: $H_{t}$' )
    xlabel( 'Normal Load: $F_{z}$' )

    linkaxes( [ax1, ax2, ax3, ax4], 'x' )
    
    subplot(20,2,33:2:39)
    fsurf(@(Fz, Slip) t0(12, Fz, 0, Slip), [0 2500 -15 15])
    ylabel( 'Shift: $H_{t}$' )
    xlabel( 'Normal Load: $F_{z}$' )
    zlabel( 'Trail: $t_{0}$' )
    
    ax5 = subplot(20,2,2:2:10);
    fplot( @(Fz) Br(10, Fz, 0), [0 2500], 'b'   ); hold on;
    fplot( @(Fz) Br(12, Fz, 0), [0 2500], 'm'   );
    fplot( @(Fz) Br(12, Fz, 2), [0 2500], 'm--' );
    fplot( @(Fz) Br(12, Fz, 4), [0 2500], 'm:'  );
    fplot( @(Fz) Br(14, Fz, 0), [0 2500], 'r'   );

    scatter( [Mesh(:).Load], [Nominal(:).qbz10] .* ...
        By([Mesh(:).Pressure], [Mesh(:).Load], [Mesh(:).Inclination]), 'k.' ); 
    ylabel( 'Curvature: $B_{r}$' )
    xlabel( 'Normal Load: $F_{z}$' )
    
    legend( {'$P = 10$ [$psi$], $\gamma = 0$ [$deg$]', ...
             '$P = 12$ [$psi$], $\gamma = 0$ [$deg$]', ...
             '$P = 12$ [$psi$], $\gamma = 2$ [$deg$]', ...
             '$P = 12$ [$psi$], $\gamma = 4$ [$deg$]', ...
             '$P = 14$ [$psi$], $\gamma = 0$ [$deg$]', ...
             'Nominal Fit Coefficient' }, 'Interpreter', 'latex' );
         
    ax6 = subplot(20,2,12:2:20);
    fplot( @(Fz) Dr(10, Fz, 0), [0 2500], 'b'   ); hold on;
    fplot( @(Fz) Dr(12, Fz, 0), [0 2500], 'm'   );
    fplot( @(Fz) Dr(12, Fz, 2), [0 2500], 'm--' );
    fplot( @(Fz) Dr(12, Fz, 4), [0 2500], 'm:'  );
    fplot( @(Fz) Dr(14, Fz, 0), [0 2500], 'r'   );

    scatter( [Mesh(:).Load], [Nominal(:).Dr], 'k.' ); 
    ylabel( 'Peak: $D_{r}$' )
    xlabel( 'Normal Load: $F_{z}$' )          

    ax7 = subplot(20,2,22:2:30);
    fplot( @(Fz) Hf(10, Fz, 0), [0 2500], 'b'   ); hold on;
    fplot( @(Fz) Hf(12, Fz, 0), [0 2500], 'm'   );
    fplot( @(Fz) Hf(12, Fz, 2), [0 2500], 'm--' );
    fplot( @(Fz) Hf(12, Fz, 4), [0 2500], 'm:'  );
    fplot( @(Fz) Hf(14, Fz, 0), [0 2500], 'r'   );

    scatter( [Mesh(:).Load], [Nominal(:).Hf], 'k.' ); 
    ylabel( 'Shift: $H_{f}$' )
    xlabel( 'Normal Load: $F_{z}$' )    
    
    linkaxes( [ax5, ax6, ax7], 'x' );
    
    subplot(20,2,32:2:40)
    fsurf(@(Fz, Slip) Mzro(12, Fz, 0, Slip), [0 2500 -15 15])
    ylabel( 'Shift: $H_{t}$' )
    xlabel( 'Normal Load: $F_{z}$' )
    zlabel( 'Residual: $M_{zr0}$' )

    sgtitle( 'Pure Aligning Moment ($M_{z0}$) Response Surfaces' )
end
        
%% Surface Plotting
Figure.Mzo.Surfaces = figure( 'Name'       , 'Aligning Moment Surfaces', ...
                              'NumberTitle', 'off', ...
                              'Visible'    , 'on' );
    
for p = 1 : size( Raw, 1 )
    for c = 1 : size( Raw, 3 )
        subplot( size( Raw, 3 ), size( Raw, 1 ), ...
            sub2ind( [size( Raw, 1 ), size( Raw, 3 )], p, c ) );
        
        plot3( [Raw(p,:,c).Load], rad2deg([Raw(p,:,c).Slip]), [Raw(p,:,c).Moment], 'k.' ); hold on;
        fsurf( @(Fz, Slip) Mzo( Mesh(p,1,c).Pressure, Fz, ...
            Mesh(p,1,c).Inclination, Slip ), [0 2500 -15 15] )
        
        xlabel( 'Normal Load ($F_{z}$) [$N$]' )
        ylabel( 'Slip Angle ($\alpha$) [$deg$]' )
        zlabel( 'Aligning Moment ($M_{z}$) [$N-m$]' )
        title( { ['Pressure ($P_{i}$): $', num2str(Mesh(p,1,c).Pressure), '$ [$psi$]'], ...
                ['Inclination ($\gamma$): $', num2str(Mesh(p,1,c).Inclination), '$ [$deg$]'] } )
    end
end

sgtitle( 'Pure Aligning MF6.1 Pacejka Fit' )
Figure.Mzo.Surfaces.WindowState = 'minimized';

%% Local Functions
    function [By, Bt, Dt, Et, Ht, Br, Dr, Hf, t0, Mzro, Mzo] = VariantEval( Tire )
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
        
        % Pneumatic Trail & Residual Moment
        Bt = @(Fz, Gam) (Tire.Pacejka.q.B.z(1) + Tire.Pacejka.q.B.z(2).*dFz(Fz) + Tire.Pacejka.q.B.z(3).*dFz(Fz).^2) .* ...
            (1 + Tire.Pacejka.q.B.z(5).*abs(Gam) + Tire.Pacejka.q.B.z(6).*Gam.^2);
        
        Ct = Tire.Pacejka.q.C.z(1);
        
        Dt = @(Pi, Fz, Gam) Tire.Pacejka.Ro .* (Fz./Tire.Pacejka.Fzo) .* ...
            ( Tire.Pacejka.q.D.z(1) + Tire.Pacejka.q.D.z(2).*dFz(Fz) ) .* ( 1 - Tire.Pacejka.p.P.z(1) .*dPi(Pi) ) .* ...
            ( 1 + Tire.Pacejka.q.D.z(3).*abs(Gam) + Tire.Pacejka.q.D.z(4).*Gam.^2 );
        
        Ht = @(Fz, Gam) Tire.Pacejka.q.H.z(1) + Tire.Pacejka.q.H.z(2).*dFz(Fz) + ...
            (Tire.Pacejka.q.H.z(3) + Tire.Pacejka.q.H.z(4).*dFz(Fz)).*Gam;
        
        Et = @(Fz, Gam, Slip) (Tire.Pacejka.q.E.z(1) + Tire.Pacejka.q.E.z(2).*dFz(Fz) + Tire.Pacejka.q.E.z(3).*dFz(Fz).^2) .* ...
            (1 + (Tire.Pacejka.q.E.z(4) + Tire.Pacejka.q.E.z(5).*Gam).*(2/pi).* ...
            atan( Bt(Fz, Gam) .* Ct .* deg2rad(Slip) + Ht(Fz, Gam) ) );
        
        Br = @(Pi, Fz, Gam) Tire.Pacejka.q.B.z(10) .* By(Pi, Fz, Gam) .* Cy;
        
        Cr = 1;
        
        Dr = @(Pi, Fz, Gam) Tire.Pacejka.Ro .* Fz .* ( ( Tire.Pacejka.q.D.z(6) + Tire.Pacejka.q.D.z(7).*dFz(Fz) ) + ...
            ( ( Tire.Pacejka.q.D.z(8) + Tire.Pacejka.q.D.z(9).*dFz(Fz) ) .* ( 1 + Tire.Pacejka.p.P.z(2).*dPi(Pi) ) + ...
            ( Tire.Pacejka.q.D.z(10) + Tire.Pacejka.q.D.z(11).*dFz(Fz) ) .* abs(Gam) ).*Gam );
        
        Hf = @(Pi, Fz, Gam) Hy(Pi, Fz, Gam) + Vy(Fz, Gam) ./ Kya(Pi, Fz, Gam);
        
        % Evaluate Functions & Error
        t0 = @(Pi, Fz, Gam, Slip) Dt(Pi, Fz, Gam) .* cos( Ct .* ...
            atan( (1-Et(Fz, Gam, Slip)).*(Bt(Fz, Gam).*(deg2rad(Slip)+Ht(Fz, Gam)) + ...
            Et(Fz, Gam, Slip).*atan( Bt(Fz, Gam).*(deg2rad(Slip)+Ht(Fz,Gam)) ) ) ) ) .* cos( deg2rad(Slip) );
        
        Mzro = @(Pi, Fz, Gam, Slip) Dr(Pi, Fz, Gam) .* cos( Cr .* ...
            atan( Br(Pi, Fz, Gam).*(deg2rad(Slip)+Hf(Pi, Fz, Gam)) ) ) .* cos( deg2rad(Slip) );
        
        Mzo = @(Pi, Fz, Gam, Slip) -t0(Pi, Fz, Gam, Slip) .* ...
            Fyo(Pi, Fz, Gam, Slip) + Mzro(Pi, Fz, Gam, Slip);
    end
end