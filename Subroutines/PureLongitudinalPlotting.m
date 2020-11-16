function [] = PureLongitudinalPlotting( Mesh, Raw, Nominal, Response, Variant, Tire )
% Plots Results from Pure Longitudinal (Fxo) MF6.1 Fitting

%% Declare Global Variables
global Figure

%% Nominal Plotting
for p = 1 : size( Nominal, 1 )
    Figure.Fxo.Nominal(p) = figure( ...
        'Name', ['Nominal Fits, Pressure (p): ', num2str(Mesh(p,1,1).Pressure), ' [psi]'], ...
        'NumberTitle', 'off');
    
    for z = 1 : size( Nominal, 2 )
        for c = 1 : size( Nominal, 3 )
            subplot( size( Nominal, 3 ), size( Nominal, 2 ), ...
                sub2ind( [size( Nominal, 2 ), size( Nominal, 3 )], z, c ) );
        
            plot( Raw(p,z,c).Slip, Raw(p,z,c).Force, 'k.' ); hold on;
            plot( Raw(p,z,c).Slip, Nominal(p,z,c).Residual, 'r.')
            fplot( NominalEval( Nominal(p,z,c).C, Nominal(p,z,c).D, ...
                Nominal(p,z,c).E, Nominal(p,z,c).K, ...
                Nominal(p,z,c).H, Nominal(p,z,c).V ), ...
                [-0.25 0.25], 'g' )
            
            xlabel( 'Slip Ratio ($\kappa$) [$deg$]' )
            ylabel( 'Longitudinal Force ($F_{x}$) [$N$]' )
            title( { ['Normal Load ($F_{z}$): $', num2str(round(Mesh(p,z,c).Load,1)), '$ [$N$]'], ...
                ['Camber ($\gamma$): $', num2str(Mesh(p,z,c).Camber), '$ [$deg$]'] } )
        end
    end
    
    sgtitle( {'Nominal P6 Pacejka Fits', ...
        ['Pressure ($P_{i}$): $', num2str(Mesh(p,1,1).Pressure), '$ [$psi$]'] } )
    
    Figure.Fxo.Nominal(p).WindowState = 'minimized';
end

%% Response Surface Plotting
Response; %#ok<VUNUS>
Variant; %#ok<VUNUS>

%% Variant Plotting
Fxo = VariantEval( Tire );

Figure.Fxo.Variant = figure( 'Name', 'Variant Fit', 'NumberTitle', 'off');
    
for p = 1 : size( Raw, 1 )
    for c = 1 : size( Raw, 3 )
        subplot( size( Raw, 3 ), size( Raw, 1 ), ...
            sub2ind( [size( Raw, 1 ), size( Raw, 3 )], p, c ) );
        
        plot3( [Raw(p,:,c).Load], [Raw(p,:,c).Slip], [Raw(p,:,c).Force], 'k.' ); hold on;
        fsurf( @(Fz, Slip) Fxo( Mesh(p,1,c).Pressure, Fz, ...
            Mesh(p,1,c).Camber, Slip ), [0 2500 -0.25 0.25] )
        
        xlabel( 'Normal Load ($F_{z}$) [$N$]' )
        ylabel( 'Slip Ratio ($\kappa$) [$deg$]' )
        zlabel( 'Longitudinal Force ($F_{x}$) [$N$]' )
        title( { ['Pressure ($P_{i}$): $', num2str(Mesh(p,1,c).Pressure), '$ [$psi$]'], ...
                ['Camber ($\gamma$): $', num2str(Mesh(p,1,c).Camber), '$ [$deg$]'] } )
    end
end

sgtitle( 'Pure Longitudinal MF6.1 Pacejka Fit' )
Figure.Fxo.Variant.WindowState = 'minimized';

%% Local Functions
    function FxoNom = NominalEval( C, D, E, K, H, V )
        % Magic Formula Evaluation
        B = K ./ ( C.*D );
        
        FxoNom = @(Slip) ( D.*sin( C.*atan( ...
            ( 1 - E ).*( B.*( Slip - H ) ) + ...
            E.*atan( B .* ( Slip - H ) ) ) ) ) + V;
    end

    function FxoVar = VariantEval( Tire )
        dFz = @(Fz) (Fz - Tire.Fzo) ./ Tire.Fzo;
        dPi = @(Pi) (Pi - Tire.Pio) ./ Tire.Pio;
        
        Cx = Tire.p.C.x(1);
        
        Dx = @(Pi, Fz, Gam) (Tire.p.D.x(1) + Tire.p.D.x(2).*dFz(Fz)) .* ...
            (1 + Tire.p.P.x(3).*dPi(Pi) + Tire.p.P.x(4).*dPi(Pi).^2) .* ...
            (1 - Tire.p.D.x(3).*Gam.^2).*Fz;
        
        Ex = @(Fz, Slip) ( Tire.p.E.x(1) + Tire.p.E.x(2).*dFz(Fz) ...
            + Tire.p.E.x(3).*dFz(Fz).^2 ) .* ( 1 - Tire.p.E.x(4).*sign(Slip) );
        
        Kxk = @(Pi, Fz) Fz.*(Tire.p.K.x(1) + Tire.p.K.x(2).*dFz(Fz) ) .* ...
            exp( Tire.p.K.x(3) .* dFz(Fz) ) .* ...
            (1 + Tire.p.P.x(1).*dPi(Pi) + Tire.p.P.x(2).*dPi(Pi).^2);
        
        Bx = @(Pi, Fz, Gam) Kxk(Pi, Fz) ./ ( Cx.*Dx(Pi, Fz, Gam) );
        
        Vx = @(Fz) Fz.*(Tire.p.V.x(1) + Tire.p.V.x(2).*dFz(Fz) );
        
        Hx = @(Fz) Tire.p.H.x(1) + Tire.p.H.x(2).*dFz(Fz);
        
        FxoVar = @(Pi, Fz, Gam, Slip) Dx(Pi, Fz, Gam) .* ...
            sin( Cx .* atan( (1-Ex(Fz, Slip)) .* ...
            Bx(Pi, Fz, Gam).*(Slip - Hx(Fz) ) + ...
            Ex(Fz, Slip).*atan( ...
            Bx(Pi, Fz, Gam).*(Slip - Hx(Fz) ) ) ) ) + Vx(Fz);
    end
end