function OverturningPlotting( Mesh, Raw, Tire, Figure )
%% Overturning Plotting - Plots Results from Overturning (Mx) Fitting
% Plots variant values from fitting process for overturning process for the
% equations given in Pacejka's "Tire and Vehicle Dynamics" [3rd Edition] in
% section 4.3.2 (page 176). 
%
% Inputs: 
%   Mesh    - Operating Condition Space
%   Raw     - Allocated Data
%   Variant - Parameter Fit Values for Pacejka Coefficients
%   Tire    - Tire Model
%   Figure  - Stores Model Figures
%
% Outputs:
%   - Surface Plots of Overturning Fitting
%
% Authors:
%   Carlos Lopez (calopez@ucdavis.edu) [Dec 2020 - June 2022]
%
% Last Updated: 01-May-2021

%% Variant Surface Plotting
Figure.Mx.Surfaces = figure( 'Name'        , 'Overturning Moment Surfaces', ...
                              'NumberTitle', 'off', ...
                              'Visible'    , 'on' );

Alpha = -20:20;
Load  = linspace(0,2600,20);
[Load, Alpha] = meshgrid( Load, Alpha );

for p = 1 : size( Raw, 1 )
    for c = 1 : size( Raw, 3 )
        subplot( size( Raw, 3 ), size( Raw, 1 ), ...
            sub2ind( [size( Raw, 1 ), size( Raw, 3 )], p, c ) );
        
        plot3( [Raw(p,:,c).Load], rad2deg([Raw(p,:,c).Alpha]), [Raw(p,:,c).Moment], 'k.' ); hold on;
        surf( Load, Alpha, VariantEval( Tire, Alpha, 0, Load, Mesh(p,1,c).Pressure, ...
            Mesh(p,1,c).Inclination) )
    
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
function [Mx] = VariantEval( Tire, Alpha, Kappa, Load, Pressure, Inclination )
    [~, ~, ~, Mx, ~] = ContactPatchLoads(Tire, Alpha, Kappa, ...
        Load, Pressure, Inclination, 10, 1, ...
        struct('Pure', 'Pacejka', 'Combined', 'MNC'));
end

end