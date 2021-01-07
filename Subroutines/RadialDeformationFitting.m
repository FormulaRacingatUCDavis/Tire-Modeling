function Tire = RadialDeformationFitting( Tire, Data )
%% Parse Data
RlLoad = Data

%% Loaded Radius Fit
[xData, yData, zData] = prepareSurfaceData( RlLoad, RlPressure, RlRadius );

% Set up fittype and options.
ft = fittype( 'poly23' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Lower = [-Inf -Inf -Inf -Inf -Inf -Inf 0 -Inf 0];
opts.Upper = [ Inf  Inf  Inf  Inf  Inf  Inf 0  Inf 0];
opts.Robust = 'Bisquare';

Tire.Rl = fit( [xData, yData], zData, ft, opts );

% Plot fit with data.
subplot(1,2,1)
plot( Tire.Rl, [xData, yData], zData );
legend( {'Polynomial Fit', 'Raw Data'}, 'Location', 'NorthEast', 'Interpreter', 'latex' );
xlabel( 'Normal Load ($F_{z}$) [$N$]', 'Interpreter', 'latex' );
ylabel( 'Pressure ($P$) [$psi$]', 'Interpreter', 'latex' );
zlabel( 'Loaded Radius ($r_{l}$) [$cm$]', 'Interpreter', 'latex' );
grid on
view( 44.6, 8.2 );

%% Effective Radius
[xData, yData, zData] = prepareSurfaceData( ReLoad, ReKappa, ReRadius );

% Set up fittype and options.
ft = fittype( 'poly12' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Robust = 'Bisquare';

% Fit model to data.
Tire.Re = fit( [xData, yData], zData, ft, opts );

% Plot fit with data.
subplot(1,2,2)
plot( Re, [xData, yData], zData );
legend( {'Polynomial Fit', 'Raw Data'}, 'Location', 'NorthEast', 'Interpreter', 'latex' );
xlabel( 'Normal Load ($F_{z}$) [$N$}', 'Interpreter', 'latex' );
ylabel( 'Slip Ratio ($\kappa$) [ ]', 'Interpreter', 'latex' );
zlabel( 'Effective Radius ($r_{e}$) [$cm$]', 'Interpreter', 'latex' );
grid on
view( -83.0, 38.3 );


