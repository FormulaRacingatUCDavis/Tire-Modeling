function [] = PureLateralPlotting( Mesh, Raw, Nominal, ~, ~, Tire )
% Plots Results from Pure Lateral (Fyo) MF6.1 Fitting

%% Declare Global Variables
global Figure

%% Evaluate Variant Surface
[Dy, Ey, Kya, Kyg0, Vy, Hy, Fyo] = VariantEval( Tire );

%% Nominal Plotting
for p = 1 : size( Nominal, 1 )
    Figure.Fyo.Nominal(p) = figure( ...
        'Name'       , ['Lateral Nominal Fits, Pressure (p): ', num2str(Mesh(p,1,1).Pressure), ' [psi]'], ...
        'NumberTitle', 'off' , ...
        'Visible'    , 'on'  );
    
    for z = 1 : size( Nominal, 2 )
        for c = 1 : size( Nominal, 3 )
            subplot( size( Nominal, 3 ), size( Nominal, 2 ), ...
                sub2ind( [size( Nominal, 2 ), size( Nominal, 3 )], z, c ) );
        
            plot( Raw(p,z,c).Slip, Raw(p,z,c).Force, 'k.' ); hold on;
            
            fplot( NominalEval( Nominal(p,z,c).C, Nominal(p,z,c).D, ...
                Nominal(p,z,c).E, Nominal(p,z,c).K, ...
                Nominal(p,z,c).H, Nominal(p,z,c).V ), ...
                [-15 15], 'g-.' )
            fplot( @(Slip) Fyo( Mesh(p,z,c).Pressure, Mesh(p,z,c).Load, ...
                Mesh(p,z,c).Inclination, Slip), [-15 15], 'g' )
            
            plot( Raw(p,z,c).Slip, Nominal(p,z,c).Residual, 'r.')
            plot( Raw(p,z,c).Slip, Raw(p,z,c).Force - ...
                Fyo( Mesh(p,z,c).Pressure   , Mesh(p,z,c).Load, ...
                     Mesh(p,z,c).Inclination, Raw(p,z,c).Slip), 'y.')
            
            xlabel( 'Slip Angle ($\alpha$) [$deg$]' )
            ylabel( 'Lateral Force ($F_{y}$) [$N$]' )
            title( { ['Normal Load ($F_{z}$): $', num2str(round(Mesh(p,z,c).Load,1)), '$ [$N$]'], ...
                ['Inclination ($\gamma$): $', num2str(Mesh(p,z,c).Inclination), '$ [$deg$]'] } )
            
            if all([p z c] == ones(3,1))
                legend( {'Raw Data', 'Nominal Fit', 'Variant Fit',...
                    'Nominal Residual', 'Variant Residual'} )
            end
        end
    end
    
    sgtitle( {'Nominal P6 Pacejka Fits', ...
        ['Pressure ($P_{i}$): $', num2str(Mesh(p,1,1).Pressure), '$ [$psi$]'] } )
    
    Figure.Fyo.Nominal(p).WindowState = Figure.State;
end

%% Response Surface Plotting
if strcmpi( Figure.Mode, 'Debug' )
    Figure.Fyo.Response = figure( ...
        'Name'       , 'Lateral Response Surfaces', ...
        'NumberTitle', 'off', ...
        'Visible'    , 'on' );
    
    ax1 = subplot(6,1,1);  
    fplot( @(Fz) Dy(10, Fz, 0), [0 2500], 'b'   ); hold on;
    fplot( @(Fz) Dy(12, Fz, 0), [0 2500], 'm'   );
    fplot( @(Fz) Dy(12, Fz, 2), [0 2500], 'm--' );
    fplot( @(Fz) Dy(12, Fz, 4), [0 2500], 'm:'  );
    fplot( @(Fz) Dy(14, Fz, 0), [0 2500], 'r'   );

    scatter( [Mesh(:).Load], [Nominal(:).D], 'k.' ); 
    ylabel( 'Peak: $D_{y}$' )
    xlabel( 'Normal Load: $F_{z}$' )
    
    legend( {'$P = 10$ [$psi$], $\gamma = 0$ [$deg$]', ...
             '$P = 12$ [$psi$], $\gamma = 0$ [$deg$]', ...
             '$P = 12$ [$psi$], $\gamma = 2$ [$deg$]', ...
             '$P = 12$ [$psi$], $\gamma = 4$ [$deg$]', ...
             '$P = 14$ [$psi$], $\gamma = 0$ [$deg$]', ...
             'Nominal Fit Coefficient' }, 'Interpreter', 'latex' );
         
    ax2 = subplot(6,1,2);
    fplot( @(Fz) Ey(Fz, 0, 1, 1), [0 2500], 'm'   ); hold on;
    fplot( @(Fz) Ey(Fz, 2, 1, 1), [0 2500], 'm--' );
    fplot( @(Fz) Ey(Fz, 4, 1, 1), [0 2500], 'm:'  );
    scatter( [Mesh(:).Load], [Nominal(:).E], 'k.' ); 
    ylabel( 'Curvature: $E_{y}$' )
    xlabel( 'Normal Load: $F_{z}$' )
    
    ax3 = subplot(6,1,3);
    fplot( @(Fz) Kya(12, Fz, 0), [0 2500], 'm'   ); hold on;
    fplot( @(Fz) Kya(12, Fz, 2), [0 2500], 'm--' );
    fplot( @(Fz) Kya(12, Fz, 4), [0 2500], 'm:'  );
    fplot( @(Fz) Kya(10, Fz, 0), [0 2500], 'b'   );
    fplot( @(Fz) Kya(14, Fz, 0), [0 2500], 'r'   );
    scatter( [Mesh(:).Load], [Nominal(:).K], 'k.' ); 
    ylabel( 'Stiffness: $K_{y \alpha}$' )
    xlabel( 'Normal Load: $F_{z}$' )
    
    ax4 = subplot(6,1,4);
    fplot( @(Fz) Hy(12, Fz, 0), [0 2500], 'm'   ); hold on;
    fplot( @(Fz) Hy(12, Fz, 2), [0 2500], 'm--' );
    fplot( @(Fz) Hy(12, Fz, 4), [0 2500], 'm:'  );
    scatter( [Mesh(:).Load], [Nominal(:).H], 'k.' ); 
    ylabel( 'Horizontal: $H_{y}$' )
    xlabel( 'Normal Load: $F_{z}$' )
    
    ax5 = subplot(6,1,5);
    fplot( @(Fz) Vy(Fz,0), [0 2500], 'm'   ); hold on;
    fplot( @(Fz) Vy(Fz,0), [0 2500], 'm--' );
    fplot( @(Fz) Vy(Fz,0), [0 2500], 'm:'  );
    scatter( [Mesh(:).Load], [Nominal(:).V], 'k.' ); 
    ylabel( 'Vertical: $V_{y}$' )
    xlabel( 'Normal Load: $F_{z}$' )
    
    ax6 = subplot(6,1,6);
    fplot( @(Fz) Kyg0(12,Fz), [0 2500], 'm' ); hold on;
    fplot( @(Fz) Kyg0(10,Fz), [0 2500], 'b' );
    fplot( @(Fz) Kyg0(14,Fz), [0 2500], 'r' );
    ylabel( 'Camber: $K_{y \gamma 0}$' )
    xlabel( 'Normal Load: $F_{z}$' )
    
    linkaxes( [ax1, ax2, ax3, ax4, ax5, ax6], 'x' )
        
    Figure.Fyo.Response.WindowState = Figure.State;
end

%% Variant Surface Plotting
Figure.Fyo.Surfaces = figure( 'Name', 'Lateral Force Surfaces', 'NumberTitle', 'off', 'Visible', 'on' );

for p = 1 : size( Raw, 1 )
    for c = 1 : size( Raw, 3 )
        subplot( size( Raw, 3 ), size( Raw, 1 ), ...
            sub2ind( [size( Raw, 1 ), size( Raw, 3 )], p, c ) );
        
        plot3( [Raw(p,:,c).Load], [Raw(p,:,c).Slip], [Raw(p,:,c).Force], 'k.' ); hold on;
        fsurf( @(Fz, Slip) Fyo( Mesh(p,1,c).Pressure, Fz, ...
            Mesh(p,1,c).Inclination, Slip ), [0 2500 -15 15] )
        
        xlabel( 'Normal Load ($F_{z}$) [$N$]' )
        ylabel( 'Slip Angle ($\alpha$) [$deg$]' )
        zlabel( 'Lateral Force ($F_{y}$) [$N$]' )
        title( { ['Pressure ($P_{i}$): $', num2str(Mesh(p,1,c).Pressure), '$ [$psi$]'], ...
                ['Inclination ($\gamma$): $', num2str(Mesh(p,1,c).Inclination), '$ [$deg$]'] } )
    end
end

sgtitle( 'Pure Lateral MF6.1 Pacejka Fit' )
Figure.Fyo.Surfaces.WindowState = Figure.State;

%% Local Functions
    function Fyo = NominalEval( C, D, E, K, H, V )
        % Magic Formula Evaluation
        B = K ./ ( C.*D + 0.1);
        
        Fyo = @(Slip) ( D.*sin( C.*atan( ...
            ( 1 - E ).*( B.*( Slip + H ) ) + ...
            E.*atan( B .* ( Slip + H ) ) ) ) ) + V;
    end

    function [Dy, Ey, Kya, Kyg0, Vy, Hy, Fyo] = VariantEval( Tire )
        dFz = @(Fz) (Fz - Tire.Pacejka.Fzo) ./ Tire.Pacejka.Fzo;
        dPi = @(Pi) (Pi - Tire.Pacejka.Pio) ./ Tire.Pacejka.Pio;
        
        Cy = Tire.p.C.y(1);
        
        Dy = @(Pi, Fz, Gam) (Tire.p.D.y(1) + Tire.p.D.y(2).*dFz(Fz)) .* ...
            (1 + Tire.p.P.y(3).*dPi(Pi) + Tire.p.P.y(4).*dPi(Pi).^2) .* ...
            (1 - Tire.p.D.y(3).*Gam.^2).*Fz;
        
        Kya = @(Pi, Fz, Gam) Tire.p.K.y(1) .* Tire.Pacejka.Fzo .* ( 1 + Tire.p.P.y(1).*dPi(Pi) ) .* ...
            ( 1 - Tire.p.K.y(3).*abs(Gam) ) .* sin( Tire.p.K.y(4) .* ...
            atan( (Fz./Tire.Pacejka.Fzo) ./ ...
            ( ( Tire.p.K.y(2) + Tire.p.K.y(5).*Gam.^2 ) .* ( 1 + Tire.p.P.y(2).*dPi(Pi) ) ) ) );
        
        Kyg0 = @(Pi, Fz) Fz.*(Tire.p.K.y(6) + Tire.p.K.y(7).*dFz(Fz) ) .* (1 + Tire.p.P.y(5).*dPi(Pi) );
        
        By = @(Pi, Fz, Gam) Kya(Pi, Fz, Gam) ./ ( Cy.*Dy(Pi, Fz, Gam) );
        
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
    end
end