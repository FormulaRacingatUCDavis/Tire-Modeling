function Tire = RadialDeflection( Tire, Data )
%%% Radial Deflection Modeling - Loaded & Effective Radius Model
%
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
%   Tire - FRUCDTire Object with Update Radius Property
%
% Author(s):
% Blake Christierson (bechristierson@ucdavis.edu)
%
% Last Updated: 8-Feb-2021

%% Loaded Radius
Loaded.Load     = [Data(1).Force(3,:)   , Data(2).Force(3,:)   , Data(3).Force(3,:)   ];
Loaded.Pressure = [Data(1).Pressure     , Data(2).Pressure     , Data(3).Pressure     ]; 
Loaded.Radius   = [Data(1).Radius.Loaded, Data(2).Radius.Loaded, Data(3).Radius.Loaded];

Tire.Radius.Loaded = fit( [Loaded.Load; Loaded.Pressure]', Loaded.Radius', 'poly22' );

%% Effective Radius
Form = fittype( 'sf*(p00+p10*x+p01*y+p20*x^2+p02*y^2+p11*x*y)', ...
    'independent', {'x', 'y'}, 'dependent', 'z' );

Opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
Opts.Display    = 'Off';
Opts.Lower      = [-Inf -Inf -Inf -Inf -Inf -Inf 1];
Opts.StartPoint = [ 0    0    0    0    0    0   1];
Opts.Upper      = [ Inf  Inf  Inf  Inf  Inf  Inf 1];

Effective.Load     = [Data(4).Force(3,:)      , Data(5).Force(3,:)      , Data(6).Force(3,:)      ];
Effective.Pressure = [Data(4).Pressure        , Data(5).Pressure        , Data(6).Pressure        ];
Effective.Radius   = [Data(4).Radius.Effective, Data(5).Radius.Effective, Data(6).Radius.Effective];

Tire.Radius.Effective = fit( [Effective.Load; Effective.Pressure]', Effective.Radius', Form, Opts );

%% Scaling Effective Radius
if ~all( strcmpi({Data(2:end).Tire}, Data(1).Tire ) )
    Loaded2.Load     = [Data(4).Force(3,:)   , Data(5).Force(3,:)   , Data(6).Force(3,:)   ];
    Loaded2.Pressure = [Data(4).Pressure     , Data(5).Pressure     , Data(6).Pressure     ]; 
    Loaded2.Radius   = [Data(4).Radius.Loaded, Data(5).Radius.Loaded, Data(6).Radius.Loaded];
    
    Loaded2 = fit( [Loaded2.Load; Loaded2.Pressure]', Loaded2.Radius', 'poly22' );
    
    Tire.Radius.Effective.sf = ...
        Tire.Radius.Loaded( Tire.Pacejka.Fzo./4, Tire.Pacejka.Pio ) ./ ...
        Loaded2( Tire.Pacejka.Fzo./4, Tire.Pacejka.Pio );
end

end