function PureLongitudinalPlotting( Tire, Raw, Mesh, Nominal, ~, ~, Figure )
% Plots Results from Pure Longitudinal (Fxo) MF6.1 Fitting

%% Evaluate Variant
[Dx, Ex, Kxk, Vx, Hx, Fxo] = VariantEval( Tire );

Pressure    = unique( [Mesh.Pressure]    );
Inclination = unique( [Mesh.Inclination] );

%% Response Surface Plotting
if strcmpi( Figure.Mode, 'Debug' )
    Figure.Fxo.Response = figure( ...
        'Name'       , 'Longitudinal Response Surfaces', ...
        'NumberTitle', 'off', ...
        'Visible'    , 'on' );
    
    ax1 = subplot(5,1,1); 
    fplot( @(Fz) Dx(10, Fz, 0), [0 2500], 'b'   ); hold on;
    fplot( @(Fz) Dx(12, Fz, 0), [0 2500], 'm'   );
    fplot( @(Fz) Dx(12, Fz, 2), [0 2500], 'm--' );
    fplot( @(Fz) Dx(12, Fz, 4), [0 2500], 'm:'  );
    fplot( @(Fz) Dx(14, Fz, 0), [0 2500], 'r'   );

    scatter( [Mesh(:).Load], [Nominal(:).D], 'k.' ); 
    ylabel( 'Peak: $D_{x}$' )
    xlabel( 'Normal Load: $F_{z}$' )
    
    legend( {'$P = 10$ [$psi$], $\gamma = 0$ [$deg$]', ...
             '$P = 12$ [$psi$], $\gamma = 0$ [$deg$]', ...
             '$P = 12$ [$psi$], $\gamma = 2$ [$deg$]', ...
             '$P = 12$ [$psi$], $\gamma = 4$ [$deg$]', ...
             '$P = 14$ [$psi$], $\gamma = 0$ [$deg$]', ...
             'Nominal Fit Coefficient' }, 'Interpreter', 'latex' );
         
    ax2 = subplot(5,1,2);
    fplot( @(Fz) Ex(Fz,  1), [0 2500], 'm' ); hold on;
    fplot( @(Fz) Ex(Fz, -1), [0 2500], 'm' ); hold on;
    scatter( [Mesh(:).Load], [Nominal(:).E], 'k.' ); 
    ylabel( 'Curvature: $E_{x}$' )
    xlabel( 'Normal Load: $F_{z}$' )
    
    ax3 = subplot(5,1,3);
    fplot( @(Fz) Kxk(12, Fz), [0 2500], 'm' ); hold on;
    fplot( @(Fz) Kxk(10, Fz), [0 2500], 'b' );
    fplot( @(Fz) Kxk(14, Fz), [0 2500], 'r' );
    scatter( [Mesh(:).Load], [Nominal(:).K], 'k.' ); 
    ylabel( 'Stiffness: $K_{x \kappa}$' )
    xlabel( 'Normal Load: $F_{z}$' )
    
    ax4 = subplot(5,1,4);
    fplot( @(Fz) Hx(Fz), [0 2500], 'm' ); hold on;
    scatter( [Mesh(:).Load], [Nominal(:).H], 'k.' ); 
    ylabel( 'Horizontal: $H_{x}$' )
    xlabel( 'Normal Load: $F_{z}$' )
    
    ax5 = subplot(5,1,5);
    fplot( @(Fz) Vx(Fz), [0 2500], 'm' ); hold on;
    scatter( [Mesh(:).Load], [Nominal(:).V], 'k.' ); 
    ylabel( 'Vertical: $V_{x}$' )
    xlabel( 'Normal Load: $F_{z}$' )
    
    linkaxes( [ax1, ax2, ax3, ax4, ax5], 'x' )
    
    Figure.Fxo.Response.WindowState = Figure.State;
end

%% Carpet Plots
Figure.Fxo.Surfaces = figure( 'Name'       , 'Longitudinal Force Surfaces', ...
                              'NumberTitle', 'off', ...
                              'Visible'    , 'on' );

for i = 1 : numel( Nominal )
    p = find( Mesh(i).Pressure    == Pressure    );
    c = find( Mesh(i).Inclination == Inclination );
   
    subplot( numel(Inclination), numel(Pressure), ...
        sub2ind( [numel(Pressure), numel(Inclination)], p, c ) );
    
    plot3( [Raw(i).Load], [Raw(i).Slip], [Raw(i).Force], 'k.' ); hold on;
    fsurf( @(Fz, Slip) Fxo( Mesh(i).Pressure, Fz, ...
        Mesh(i).Inclination, Slip ), [0 2500 -0.25 0.25] )
        
    xlabel( 'Normal Load ($F_{z}$) [$N$]' )
    ylabel( 'Slip Ratio ($\kappa$) []' )
    zlabel( 'Longitudinal Force ($F_{x}$) [$N$]' )
    title( { ['Pressure ($P_{i}$): $'    , num2str(Mesh(i).Pressure)   , '$ [$psi$]'], ...
             ['Inclination ($\gamma$): $', num2str(Mesh(i).Inclination), '$ [$deg$]'] } )
end

sgtitle( 'Pure Longitudinal MF6.1 Pacejka Fit' )
Figure.Fxo.Surfaces.WindowState = Figure.State;

%% Local Functions
function Fxo = NominalEval( C, D, E, K, H, V )
    % Magic Formula Evaluation
    B = K ./ ( C.*D );

    Fxo = @(Slip) ( D.*sin( C.*atan( ...
        ( 1 - E ).*( B.*( Slip + H ) ) + ...
        E.*atan( B .* ( Slip + H ) ) ) ) ) + V;
end

function [Dx, Ex, Kxk, Vx, Hx, Fxo] = VariantEval( Tire )
    dFz = @(Fz) (Fz - Tire.Pacejka.Fzo) ./ Tire.Pacejka.Fzo;
    dPi = @(Pi) (Pi - Tire.Pacejka.Pio) ./ Tire.Pacejka.Pio;

    Cx = Tire.Pacejka.p.C.x(1);

    Dx = @(Pi, Fz, Gam) (Tire.Pacejka.p.D.x(1) + Tire.Pacejka.p.D.x(2).*dFz(Fz)) .* ...
        (1 + Tire.Pacejka.p.P.x(3).*dPi(Pi) + Tire.Pacejka.p.P.x(4).*dPi(Pi).^2) .* ...
        (1 - Tire.Pacejka.p.D.x(3).*Gam.^2).*Fz;

    Ex = @(Fz, Slip) ( Tire.Pacejka.p.E.x(1) + Tire.Pacejka.p.E.x(2).*dFz(Fz) ...
        + Tire.Pacejka.p.E.x(3).*dFz(Fz).^2 ) .* ( 1 - Tire.Pacejka.p.E.x(4).*sign(Slip) );

    Kxk = @(Pi, Fz) Fz.*(Tire.Pacejka.p.K.x(1) + Tire.Pacejka.p.K.x(2).*dFz(Fz) ) .* ...
        exp( Tire.Pacejka.p.K.x(3) .* dFz(Fz) ) .* ...
        (1 + Tire.Pacejka.p.P.x(1).*dPi(Pi) + Tire.Pacejka.p.P.x(2).*dPi(Pi).^2);

    Bx = @(Pi, Fz, Gam) Kxk(Pi, Fz) ./ ( Cx.*Dx(Pi, Fz, Gam) );

    Vx = @(Fz) Fz.*(Tire.Pacejka.p.V.x(1) + Tire.Pacejka.p.V.x(2).*dFz(Fz) );

    Hx = @(Fz) Tire.Pacejka.p.H.x(1) + Tire.Pacejka.p.H.x(2).*dFz(Fz);

    Fxo = @(Pi, Fz, Gam, Slip) Dx(Pi, Fz, Gam) .* ...
        sin( Cx .* atan( (1-Ex(Fz, Slip)) .* ...
        Bx(Pi, Fz, Gam).*(Slip + Hx(Fz) ) + ...
        Ex(Fz, Slip).*atan( ...
        Bx(Pi, Fz, Gam).*(Slip + Hx(Fz) ) ) ) ) + Vx(Fz);
end

end