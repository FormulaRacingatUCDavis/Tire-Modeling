function [] = OverturningPlotting(Mesh, Raw, ~, ~, Tire)
% Plots Results from Overturning (Mx) Fitting

%% Devlare Global Variables 
global Figure

%% Evaluate Variant Surface
[Mx] = VariantEval( Tire );

%% Variant Surface Plotting
Figure.Mx.Surfaces = figure( 'Name'       , 'Lateral Force Surfaces', ...
                              'NumberTitle', 'off', ...
                              'Visible'    , 'on' );

for p = 1 : size( Raw, 1 )
    for c = 1 : size( Raw, 3 )
        subplot( size( Raw, 3 ), size( Raw, 1 ), ...
            sub2ind( [size( Raw, 1 ), size( Raw, 3 )], p, c ) );
        
        plot3( [Raw(p,:,c).Load], rad2deg([Raw(p,:,c).Slip]), [Raw(p,:,c).Force], 'k.' ); hold on;
        fsurf( @(Fz, Slip) Mx( Mesh(p,1,c).Pressure, Fz, ...
            Mesh(p,1,c).Inclination, Slip ), [0 2500 -15 15] )
        
        xlabel( 'Normal Load ($F_{z}$) [$N$]' )
        ylabel( 'Slip Angle ($\alpha$) [$deg$]' )
        zlabel( 'Lateral Force ($F_{y}$) [$N$]' )
        title( { ['Pressure ($P_{i}$): $'    , num2str(Mesh(p,1,c).Pressure)   , '$ [$psi$]'], ...
                 ['Inclination ($\gamma$): $', num2str(Mesh(p,1,c).Inclination), '$ [$deg$]'] } )
    end
end

sgtitle( 'Overturning MF6.1 Pacejka Fit' )
Figure.Mx.Surfaces.WindowState = Figure.State;

%% Local Functions

function [Mx] = VariantEval( Tire )
        dPi = @(Pi) (Pi - Tire.Pacejka.Pio) ./ Tire.Pacejka.Pio;
        
        Mx.Surface = @(Fz, Gam, Pi, x0) (Tire.Pacejka.Ro * Tire.Pacejka.Fzo) * ... 
            ( x0.qsx1 - ((x0.qsx2 * Gam) * (1 + x0.ppMx1 * dPi(Pi))) +(x0.qsx3 * ...
            (Fy(Fz, Gam, Pi, x0)./Tire.Pacejka.Fzo)) + (x0.qsx4 * cos(x0.qsx5 * ...
            atan(x0.qsx6 * (Fz./Tire.Pacejka.Fzo).^2))* sin((x0.qsx7 * Gam) + ...
            (x0.qsx8 * atan( x0.qsx9 * (Fy(Fz,Gam,Pi,x0)./Tire.Pacejka.Fzo))))) + ... 
            (x0.qsx10 * atan(x0.qsx11 * (Fz./Tire.Pacejka.Fzo)) * Gam));
end
end