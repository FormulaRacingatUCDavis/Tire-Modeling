function [] = PureLateralPlotting( Tire, Raw, Mesh, Nominal, ~, ~, Figure )
% Plots Results from Pure Lateral (Fyo) MF6.1 Fitting

%% Evaluate Variant Surface
[Dy, Ey, Kya, Kyg0, Vy, Hy, Fyo] = VariantEval( Tire );

Pressure    = unique( [Mesh.Pressure]    );
Load        = unique( [Mesh.Load]        );
Inclination = unique( [Mesh.Inclination] );

%% Nominal Plotting
for p = 1 : numel( Pressure )
    Figure.Fyo.Nominal(p) = figure( ...
        'Name'       , ['Lateral Nominal Fits, P=', num2str(Pressure(p)), ' [psi]'], ...
        'NumberTitle', 'off' , ...
        'Visible'    , 'on'  );
end

for i = 1 : numel( Nominal )
    p = find( Mesh(i).Pressure    == Pressure    );
    z = find( Mesh(i).Load        == Load        );
    c = find( Mesh(i).Inclination == Inclination );
   
    figure( Figure.Fyo.Nominal(p) )
   
    subplot( numel(Inclination), numel(Load), ...
        sub2ind( [numel(Load), numel(Inclination)], z, c ) );
    
    plot( rad2deg(Raw(i).Slip), Raw(i).Force, 'k.' ); hold on;
    
    fplot( NominalEval( Nominal(i).C, Nominal(i).D, Nominal(i).E, ...
        Nominal(i).K, Nominal(i).H, Nominal(i).V ), [-15 15], 'g-.' )
    fplot( @(Slip) Fyo( Mesh(i).Pressure, Mesh(i).Load, ...
        Mesh(i).Inclination, deg2rad(Slip) ), [-15 15], 'g' )

    plot( rad2deg(Raw(i).Slip), Nominal(i).Residual, 'r.')
    plot( rad2deg(Raw(i).Slip), Raw(i).Force - ...
        Fyo( Mesh(i).Pressure   , Mesh(i).Load, ...
             Mesh(i).Inclination, Raw(i).Slip), 'y.')

    xlabel( 'Slip Angle ($\alpha$) [$deg$]' )
    ylabel( 'Lateral Force ($F_{y}$) [$N$]' )
    title( { ['Normal Load ($F_{z}$): $', num2str(round(Mesh(i).Load,1)), '$ [$N$]'], ...
        ['Inclination ($\gamma$): $', num2str(Mesh(i).Inclination), '$ [$deg$]'] } )

    if all([z c] == ones(3,1))
        legend( {'Raw Data', 'Nominal Fit', 'Variant Fit',...
                 'Nominal Residual', 'Variant Residual'} )
    end
end

for p = 1 : numel( Pressure )
    figure( Figure.Fyo.Nominal(p) )  
    
    sgtitle( {'Nominal P6 Pacejka Fits', ...
        ['Pressure ($P_{i}$): $', num2str(Mesh(i).Pressure), '$ [$psi$]'] } )
    
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
    
    sgtitle( 'Pure Lateral Force ($F_{y0}$) Response Surfaces' )
    
    Figure.Fyo.Response.WindowState = Figure.State;
end

%% Variant Surface Plotting
Figure.Fyo.Surfaces = figure( 'Name'       , 'Lateral Force Surfaces', ...
                              'NumberTitle', 'off', ...
                              'Visible'    , 'on' );

for i = 1 : numel( Nominal )
    p = find( Mesh(i).Pressure    == Pressure    );
    c = find( Mesh(i).Inclination == Inclination );
   
    subplot( numel(Inclination), numel(Pressure), ...
        sub2ind( [numel(Pressure), numel(Inclination)], p, c ) );
    
    plot3( [Raw(i).Load], rad2deg([Raw(i).Slip]), [Raw(i).Force], 'k.' ); hold on;
    fsurf( @(Fz, Slip) Fyo( Mesh(i).Pressure, Fz, ...
        Mesh(i).Inclination, deg2rad(Slip) ), [0 2500 -15 15] )
        
    xlabel( 'Normal Load ($F_{z}$) [$N$]' )
    ylabel( 'Slip Angle ($\alpha$) [$deg$]' )
    zlabel( 'Lateral Force ($F_{y}$) [$N$]' )
    title( { ['Pressure ($P_{i}$): $'    , num2str(Mesh(i).Pressure)   , '$ [$psi$]'], ...
             ['Inclination ($\gamma$): $', num2str(Mesh(i).Inclination), '$ [$deg$]'] } )
end

sgtitle( 'Pure Lateral MF6.1 Pacejka Fit' )
Figure.Fyo.Surfaces.WindowState = Figure.State;

%% Local Functions
    function Fyo = NominalEval( C, D, E, K, H, V )
        % Magic Formula Evaluation
        B = K ./ ( C.*D );
        
        Fyo = @(Slip) ( D.*sin( C.*atan( ...
            ( 1 - E ).*( B.*( deg2rad(Slip) + H ) ) + ...
            E.*atan( B .* ( deg2rad(Slip) + H ) ) ) ) ) + V;
    end

    function [Dy, Ey, Kya, Kyg0, Vy, Hy, Fyo] = VariantEval( Tire )
        dFz = @(Fz) (Fz - Tire.Pacejka.Fzo) ./ Tire.Pacejka.Fzo;
        dPi = @(Pi) (Pi - Tire.Pacejka.Pio) ./ Tire.Pacejka.Pio;
        
        Cy = Tire.Pacejka.p.C.y(1);
        
        Dy = @(Pi, Fz, Gam) (Tire.Pacejka.p.D.y(1) + Tire.Pacejka.p.D.y(2).*dFz(Fz)) .* ...
            (1 + Tire.Pacejka.p.P.y(3).*dPi(Pi) + Tire.Pacejka.p.P.y(4).*dPi(Pi).^2) .* ...
            (1 - Tire.Pacejka.p.D.y(3).*Gam.^2).*Fz;
        
        Kya = @(Pi, Fz, Gam) Tire.Pacejka.p.K.y(1) .* Tire.Pacejka.Fzo .* ( 1 + Tire.Pacejka.p.P.y(1).*dPi(Pi) ) .* ...
            ( 1 - Tire.Pacejka.p.K.y(3).*abs(Gam) ) .* sin( Tire.Pacejka.p.K.y(4) .* ...
            atan( (Fz./Tire.Pacejka.Fzo) ./ ...
            ( ( Tire.Pacejka.p.K.y(2) + Tire.Pacejka.p.K.y(5).*Gam.^2 ) .* ( 1 + Tire.Pacejka.p.P.y(2).*dPi(Pi) ) ) ) );
        
        Kyg0 = @(Pi, Fz) Fz.*(Tire.Pacejka.p.K.y(6) + Tire.Pacejka.p.K.y(7).*dFz(Fz) ) .* ( 1 + Tire.Pacejka.p.P.y(5).*dPi(Pi) );
        
        By = @(Pi, Fz, Gam) Kya(Pi, Fz, Gam) ./ ( Cy.*Dy(Pi, Fz, Gam) );
        
        Vyg = @(Fz, Gam) Fz.*(Tire.Pacejka.p.V.y(3) + Tire.Pacejka.p.V.y(4).*dFz(Fz) ).*Gam;
        
        Vy = @(Fz, Gam) Fz.*(Tire.Pacejka.p.V.y(1) + Tire.Pacejka.p.V.y(2).*dFz(Fz) ) + Vyg(Fz, Gam);
        
        Hy = @(Pi, Fz, Gam) (Tire.Pacejka.p.H.y(1) + Tire.Pacejka.p.H.y(2).*dFz(Fz) ) .* ...
            (Kyg0(Pi, Fz).*Gam - Vyg(Fz, Gam) ) ./ Kya(Pi, Fz, Gam);
        
        Ey = @(Fz, Gam, Slip, Hy) ( Tire.Pacejka.p.E.y(1) + Tire.Pacejka.p.E.y(2).*dFz(Fz) ) .* ...
            ( 1 + Tire.Pacejka.p.E.y(5).*Gam.^2 - ...
            ( Tire.Pacejka.p.E.y(3) + Tire.Pacejka.p.E.y(4).*Gam ).*sign( Slip + Hy) );

        Fyo = @(Pi, Fz, Gam, Slip) Dy(Pi, Fz, Gam) .* ...
            sin( Cy .* atan( ( 1-Ey(Fz, Gam, Slip, Hy(Pi, Fz, Gam)) ) .* ...
            By(Pi, Fz, Gam).*( Slip + Hy(Pi, Fz, Gam) ) + ...
            Ey(Fz, Gam, Slip, Hy(Pi, Fz, Gam) ).*atan( ...
            By(Pi, Fz, Gam).*( Slip + Hy(Pi, Fz, Gam) ) ) ) ) + Vy(Fz, Gam);
    end
end