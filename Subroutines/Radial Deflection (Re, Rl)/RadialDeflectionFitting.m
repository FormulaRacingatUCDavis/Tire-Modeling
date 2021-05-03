function Tire = RadialDeflectionFitting( Tire, Data )
%% Radial Deflection Modeling - Loaded & Effective Radius Model
% This fits simple poly22 surfaces for the tire loaded and effective
% radius. Cornering data is used for determining loaded radius while drive,
% brake, & combined data is used to determine effective radius. If the
% drive, brake, & combined data is from a different tire (e.g. for 16"
% tires), the effective radius fit will be scaled.
% 
% Inputs:
%   Tire - FRUCDTire Object
%   Data - Formatted Structure of Calspan TTC Data via DataImport()
% 
% Outputs:
%   Tire - FRUCDTire Object with Updated Radius Property
%
% Author(s):
% Blake Christierson (bechristierson@ucdavis.edu)
%
% Last Updated: 02-May-2021

%% Loaded Radius
Form = fittype( 'p00+p10*Fz+p01*Pi+p20*Fz^2+p02*Pi^2+p11*Fz*Pi', ...
    'independent', {'Fz', 'Pi'}, 'dependent', 'Rl' );

Opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
Opts.Display    = 'Off';
Opts.StartPoint = [200, 0, 0, 0, 0, 0];

Loaded(1).Load     = [Data(1).Force(3,:)   , Data(2).Force(3,:)   , Data(3).Force(3,:)   ];
Loaded(1).Pressure = [Data(1).Pressure     , Data(2).Pressure     , Data(3).Pressure     ]; 
Loaded(1).Radius   = [Data(1).Radius.Loaded, Data(2).Radius.Loaded, Data(3).Radius.Loaded] .* 10;

Tire.Radius.Loaded = fit( [Loaded(1).Load; Loaded(1).Pressure]', Loaded(1).Radius', Form, Opts );

%% Effective Radius
Effective.Load      = [Data(4).Force(3,:)      , Data(5).Force(3,:)      , Data(6).Force(3,:)      ];
Effective.Pressure  = [Data(4).Pressure        , Data(5).Pressure        , Data(6).Pressure        ];
Effective.SlipRatio = [Data(4).Slip.Ratio      , Data(5).Slip.Ratio      , Data(6).Slip.Ratio      ];
Effective.Radius    = [Data(4).Radius.Effective, Data(5).Radius.Effective, Data(6).Radius.Effective] .* 10;

%%% Scaling Data Between Incoherent Data Sets
if ~all( strcmpi({Data(2:end).Tire}, Data(1).Tire ) )
    Loaded(2).Load     = [Data(4).Force(3,:)   , Data(5).Force(3,:)   , Data(6).Force(3,:)   ];
    Loaded(2).Pressure = [Data(4).Pressure     , Data(5).Pressure     , Data(6).Pressure     ]; 
    Loaded(2).Radius   = [Data(4).Radius.Loaded, Data(5).Radius.Loaded, Data(6).Radius.Loaded] .* 10;
    
    LoadedFit = fit( [Loaded(2).Load; Loaded(2).Pressure]', Loaded(2).Radius', Form, Opts );
    
    Effective.Radius = Effective.Radius .* ...
        Tire.Radius.Loaded( Tire.Pacejka.Fzo./4, Tire.Pacejka.Pio ) ./ ...
        LoadedFit( Tire.Pacejka.Fzo./4, Tire.Pacejka.Pio );
end

%%% Fit Effective Radius Surface
p000 = optimvar('p000');
p001 = optimvar('p001'); p010 = optimvar('p010'); p100 = optimvar('p100');
p011 = optimvar('p011'); p101 = optimvar('p101'); p110 = optimvar('p110');
p002 = optimvar('p002'); p020 = optimvar('p020'); p200 = optimvar('p200');

x0.p000 = 200;
x0.p001 =- 50; x0.p010 = 0; x0.p100 = 0;
x0.p011 =   0; x0.p101 = 0; x0.p110 = 0;
x0.p002 =  50; x0.p020 = 0; x0.p200 = 0;

Obj = fcn2optimexpr( @ErrorRe, p000, p001, p010, p100, p011, p101, p110, p002, p020, p200 );

[Solution, ~] = Runfmincon( Obj, x0, [], 1 );
delete( findobj( 'Type', 'figure', 'Name', 'Optimization PlotFcns' ) );

Tire.Radius.Effective = @(Kappa, Fz, Pi) Solution.p000 + ...
    Solution.p001 .* Kappa + Solution.p010 .* Fz + Solution.p100 .* Pi + ...
    Solution.p011 .* Kappa .* Fz + Solution.p101 .* Kappa .* Pi + Solution.p110 .* Fz .* Pi + ...
    Solution.p002 .* Kappa.^2 + Solution.p020 .* Fz.^2 + Solution.p200 .* Pi.^2;

%% Plotting
%%% Mesh Generation
Load      = 0:250:2000;
Pressure  = 50:10:100;
SlipRatio = linspace( -0.2, 0.2, length(Pressure) );

[~   , Pressure ] = meshgrid( Load, Pressure  );
[Load, SlipRatio] = meshgrid( Load, SlipRatio );

%%% Null Slip Ratio
NullIdx = abs(Effective.SlipRatio) < 0.005;

%%% Loaded Radius Carpet Plot (Load - Pressure)
figure( 'Name', 'Radial Deflection Models' );

subplot(2,2,[1 3])
scatter3( Loaded(1).Load, Loaded(1).Radius, Loaded(1).Pressure, [], Loaded(1).Pressure, '.', ...
    'MarkerEdgeAlpha', 0.3, 'MarkerFaceAlpha', 0.3 ); hold on;

surf( Load, Tire.Radius.Loaded( Load, Pressure ), Pressure );

xlabel( 'Normal Load, $F_{z}$ [$N$]' );
ylabel( 'Loaded Radius, $r_{l}$ [$mm$]' );
zlabel( 'Pressure, $P_{i}$ [$kPa$]' );

CB = colorbar( 'TickLabelInterpreter', 'latex' );
CB.Label.String = 'Pressure, $P_{i}$ [$kPa$]';
CB.Label.Interpreter = 'latex';

view(0,90)

%%% Effective Radius Carpet Plot (Slip Ratio - Load)
subplot(2,2,2) 
scatter3( Effective.SlipRatio, Effective.Radius, Effective.Load, [], Effective.Pressure, '.', ...
    'MarkerEdgeAlpha', 0.3, 'MarkerFaceAlpha', 0.3 ); hold on;

surf( SlipRatio, Tire.Radius.Effective( SlipRatio, Load, 70 ), Load, 70.*ones(size(Load)) );

xlabel( 'Slip Ratio, $\kappa$ [ ]' );
ylabel( 'Effective Radius, $r_{e}$ [$mm$]' );
zlabel( 'Normal Load, $F_{z}$ [$N$]' );

CB = colorbar( 'TickLabelInterpreter', 'latex' );
CB.Label.String = 'Pressure, $P_{i}$ [$kPa$]';
CB.Label.Interpreter = 'latex';

view(0,90)

%%% Effective Radius Carpet Plot (Load - Pressure)
subplot(2,2,4) 
scatter3( Effective.Load(NullIdx), Effective.Radius(NullIdx), ...
    Effective.Pressure(NullIdx), [], Effective.Pressure(NullIdx), '.', ...
    'MarkerEdgeAlpha', 0.3, 'MarkerFaceAlpha', 0.3 ); hold on;

surf( Load, Tire.Radius.Effective( 0, Load, Pressure ), Pressure );

xlabel( 'Normal Load, $F_{z}$ [$N$]' );
ylabel( 'Effective Radius, $r_{e}$ [$mm$]' );
zlabel( 'Pressure, $P_{i}$ [$kPa$]' );

CB = colorbar( 'TickLabelInterpreter', 'latex' );
CB.Label.String = 'Pressure, $P_{i}$ [$kPa$]';
CB.Label.Interpreter = 'latex';

view(0,90)

%% Local Functions
function RMSE = ErrorRe( p000, p001, p010, p100, ...
        p011, p101, p110, p002, p020, p200 )
    
    Re = p000 + p001.*Effective.SlipRatio + p010.*Effective.Load + p100.*Effective.Pressure + ...
        p011.*Effective.SlipRatio.*Effective.Load + p101.*Effective.SlipRatio.*Effective.Pressure + ...
        p110.*Effective.Load.*Effective.Pressure + ...
        p002.*Effective.SlipRatio.^2 + p020.*Effective.Load.^2 + p200.*Effective.Pressure.^2;
     
    RMSE = sqrt( mean( (Re(:) - Effective.Radius(:)).^2 ) );
end

end