function [] = PureLateralPlotting( Mesh, Raw, Nominal, Tire )
% Plots Results from Pure Lateral (Fyo) MF6.1 Fitting

%% Declare Global Variables
global Figure

%% Evaluate Variant Fit
Fyo = VariantEval( Tire );

%% Operating Case Plots
Figure.Fyo.Cases = figure(); close(gcf);

for p = 1 : size( Nominal, 1 )
    Figure.Fyo.Cases(p) = figure( ...
        'Name', ['Lateral Force Cases, Pressure (p): ', num2str(Mesh(p,1,1).Pressure), ' [psi]'], ...
        'NumberTitle', 'off');
    
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
                Mesh(p,z,c).Camber, Slip), [-15 15], 'g' )
            
            plot( Raw(p,z,c).Slip, Nominal(p,z,c).Residual, 'r.')
            plot( Raw(p,z,c).Slip, ...
                Raw(p,z,c).Force-Fyo( Mesh(p,z,c).Pressure, Mesh(p,z,c).Load, Mesh(p,z,c).Camber, Raw(p,z,c).Slip), 'y.')
            
            xlabel( 'Slip Angle ($\alpha$) [$deg$]' )
            ylabel( 'Lateral Force ($F_{y}$) [$N$]' )
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
    
    Figure.Fyo.Cases(p).WindowState = 'minimized';
end

%% Surface Plotting
Figure.Fyo.Surfaces = figure( 'Name', 'Lateral Force Surfaces', 'NumberTitle', 'off');
    
for p = 1 : size( Raw, 1 )
    for c = 1 : size( Raw, 3 )
        subplot( size( Raw, 3 ), size( Raw, 1 ), ...
            sub2ind( [size( Raw, 1 ), size( Raw, 3 )], p, c ) );
        
        plot3( [Raw(p,:,c).Load], [Raw(p,:,c).Slip], [Raw(p,:,c).Force], 'k.' ); hold on;
        fsurf( @(Fz, Slip) Fyo( Mesh(p,1,c).Pressure, Fz, ...
            Mesh(p,1,c).Camber, Slip ), [0 2500 -15 15] )
        
        xlabel( 'Normal Load ($F_{z}$) [$N$]' )
        ylabel( 'Slip Angle ($\alpha$) [$deg$]' )
        zlabel( 'Lateral Force ($F_{y}$) [$N$]' )
        title( { ['Pressure ($P_{i}$): $', num2str(Mesh(p,1,c).Pressure), '$ [$psi$]'], ...
                ['Camber ($\gamma$): $', num2str(Mesh(p,1,c).Camber), '$ [$deg$]'] } )
    end
end

sgtitle( 'Pure Lateral MF6.1 Pacejka Fit' )
Figure.Fyo.Surfaces.WindowState = 'minimized';

%% Local Functions
    function FyoNom = NominalEval( C, D, E, K, H, V )
        % Magic Formula Evaluation
        B = K ./ ( C.*D + 0.1);
        
        FyoNom = @(Slip) ( D.*sin( C.*atan( ...
            ( 1 - E ).*( B.*( Slip + H ) ) + ...
            E.*atan( B .* ( Slip + H ) ) ) ) ) + V;
    end

    function FyoVar = VariantEval( Tire )
        dFz = @(Fz) (Fz - Tire.Fzo) ./ Tire.Fzo;
        dPi = @(Pi) (Pi - Tire.Pio) ./ Tire.Pio;
        
        Cy = Tire.p.C.y(1);
        
        Dy = @(Pi, Fz, Gam) (Tire.p.D.y(1) + Tire.p.D.y(2).*dFz(Fz)) .* ...
            (1 + Tire.p.P.y(3).*dPi(Pi) + Tire.p.P.y(4).*dPi(Pi).^2) .* ...
            (1 - Tire.p.D.y(3).*Gam.^2).*Fz;
        
        Kya = @(Pi, Fz, Gam) Tire.p.K.y(1) .* Tire.Fzo .* ( 1 + Tire.p.P.y(1).*dPi(Pi) ) .* ...
            ( 1 - Tire.p.K.y(3).*abs(Gam) ) .* sin( Tire.p.K.y(4) .* ...
            atan( (Fz./Tire.Fzo) ./ ...
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

        FyoVar = @(Pi, Fz, Gam, Slip) Dy(Pi, Fz, Gam) .* ...
            sin( Cy .* atan( (1-Ey(Fz, Gam, Slip, Hy(Pi, Fz, Gam) )) .* ...
            By(Pi, Fz, Gam).*(Slip + Hy(Pi, Fz, Gam) ) + ...
            Ey(Fz, Gam, Slip, Hy(Pi, Fz, Gam) ).*atan( ...
            By(Pi, Fz, Gam).*(Slip + Hy(Pi, Fz, Gam) ) ) ) ) + Vy(Fz, Gam);
        
 %       if strcmpi( Mode, 'Debug' )
 %           Figure.Fyo.Debug = figure;
%            subplot(1,1)
             
%        end
    end
end