function [] = CombinedAligningPlotting(Mesh, Raw, Nominal, Tire)
%Plots results from Combined Lateral (MZ) MF6.1 Fitting
%% Declare Global Variables 
global Figure
%% Evaluate Variant Fit
MZ = VariantEval(Tire);

%% Operating Case Plots
Figure.MZ.Cases = figure(); close(gcf);

for p = 1 : size( Nominal, 1 )
    Figure.Fyo.Cases(p) = figure( ...
        'Name', ['Combined Aligning Cases, Pressure (p): ', num2str(Mesh(p,1,1).Pressure), ' [psi]'], ...
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
Figure.MZ.Surfaces = figure('Name', 'Combined Aligning Surfaces', 'NumberTitle', 'off');
for p = 1 : size( Raw, 1 )
    for c = 1 : size( Raw, 3 )
        subplot( size( Raw, 3 ), size( Raw, 1 ), ...
            sub2ind( [size( Raw, 1 ), size( Raw, 3 )], p, c ) );
        
        plot3( [Raw(p,:,c).Load], [Raw(p,:,c).Slip], [Raw(p,:,c).Force], 'k.' ); hold on;
        fsurf( @(Fz, Slip) MZ( Mesh(p,1,c).Pressure, Fz, ...
            Mesh(p,1,c).Camber, Slip ), [0 2500 -15 15] )
        
        xlabel( 'Normal Load ($F_{z}$) [$N$]' )
        ylabel( 'Slip Angle ($\alpha$) [$deg$]' )
        zlabel( 'Aligning Moment ($F_{y}$) [$N$]' )
        title( { ['Pressure ($P_{i}$): $', num2str(Mesh(p,1,c).Pressure), '$ [$psi$]'], ...
                ['Camber ($\gamma$): $', num2str(Mesh(p,1,c).Camber), '$ [$deg$]'] } )
    end
end

sgtitle( 'Combined Lateral MF6.1 Pacejka Fit' )
Figure.MZ.Surfaces.WindowState = 'minimized';
%% Getting Solutions from Combined Lateral and Longitudinal
Fx = abs(Fx0 .* Fy0 ./ sqrt( (Kappa - Kappa0).^2 .* Fy0.^2 + Fx0.^2 .* tan(Alpha - Alpha0).^2 ) .* ...
                        sqrt( (Kappa - Kappa0).^2 .* Kya.^2 + (1 - abs(Kappa - Kappa0)).^2 .* cos(Alpha - Alpha0).^2 .* Fx0.^2 ) ./ ...
                        Kya) .* sign(Fx0);
                    
Fy = abs( Fx0 .* Fy0 ./ sqrt( (Kappa - Kappa0).^2 .* Fy0.^2 + Fx0.^2 .* tan(Alpha - Alpha0).^2 ) .* ...
                        sqrt( (1 - abs(Kappa - Kappa0)).^2 .* cos(Alpha - Alpha0).^2 .* Fy0.^2 + sin(Alpha - Alpha0).^2 .* Kxk.^2 ) ./ ...
                        ( Kxk .* cos(Alpha - Alpha0)) ) .* sign(Fy0);

%% Local Functions 
    function MZNom = NominalEval()
        cos_prime = vcx./(vc + 0.1);
        s = r_0 * (ssz1 + (ssz2 * (Fy_0./Fz_0
        MZNom = @(Slip, Kappa) (-Dcos(C*atan(B*Slip) - E*((B*Slip) - ... 
            atan(B*Slip)))*cos_prime * Fy) + (D * cos(C*atan(B * alpha))) ...
            + 
    end
    function MZVar = VariantEval(Tire)
        dfz = x;
    end










end