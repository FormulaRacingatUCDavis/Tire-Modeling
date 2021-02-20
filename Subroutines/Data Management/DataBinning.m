function Bin = DataBinning( Data )
%% DataBinning - Separates Experimental Operating Conditions
% This function creates logical bins to separate the total data in the
% separate experimental conditions.
% 
% Inputs:
%   Data - Parsed FSAE TTC Data
% 
% Outputs:
%   Bin  - Logical Binnings for Separating Operating Conditions
%
% Author(s): 
% Blake Christierson (bechristierson@ucdavis.edu) [Sep 2018 - Jun 2021] 
% Carlos Lopez       (calopez@ucdavis.edu       ) [Jan 2019 -         ]
% Leonardo Howard    (leohoward@ucdavis.edu     ) [Feb 2021 -         ]
% 
% Last Updated: 17-Feb-2021

%%% Conversion Functions
lbf2N   = @(lbf) lbf .* 4.44822;
mph2mps = @(mph) mph .* 1.60934 ./ 3.6;
psi2kPa = @(psi) psi .* 6.89476;

%%% Bin Values
% Derived from Round 6-8 FSAE TTC Content Guides
Values.Pressure    = psi2kPa([  8,  10,  12  14,        ]'); % Pressure Bins     [kPa]
Values.Load        = lbf2N(  [ 50, 100, 150 200, 250 350]'); % Normal Force Bins [N]
Values.Inclination =         [- 2,   0,   2   4,        ]' ; % Inclination Bins  [deg]
Values.Velocity    = mph2mps([  0,   2,  15  25,  45    ]'); % Belt Speed Bins   [m/s]
Values.Slip.Angle  = deg2rad([- 1,   0,   1,  6         ]'); % Slip Angle Bins   [deg]
Values.Slip.Ratio  =            0;                           % Slip Ratio Bins   [ ]

Tolerance.Pressure    = min( diff(Values.Pressure   ) ) / 3;  % Pressure Tolerance      [psi]
Tolerance.Load        = min( diff(Values.Load       ) ) / 3;  % Normal Force Tolerance  [N]
Tolerance.Inclination = min( diff(Values.Inclination) ) / 3;  % Inclination Tolerance   [deg]
Tolerance.Velocity    = min( diff(Values.Velocity   ) ) / 3;  % Belt Velocity Tolerance [m/s]
Tolerance.Slip.Angle  = min( diff(Values.Slip.Angle ) ) / 10; % Slip Angle Tolerance    [deg]
Tolerance.Slip.Ratio  = Tolerance.Slip.Angle     * 0.15 / 12; % Slip Ratio Tolerance    [ ]

% Bin Creation
Bin.Pressure = ( Data.Pressure > (Values.Pressure - Tolerance.Pressure) ) & ...
               ( Data.Pressure < (Values.Pressure + Tolerance.Pressure) );

Bin.Load = ( Data.Force(3,:) > (Values.Load - Tolerance.Load) ) & ...
           ( Data.Force(3,:) < (Values.Load + Tolerance.Load) );

Bin.Inclination = ( Data.Inclination > Values.Inclination - Tolerance.Inclination ) & ...
                  ( Data.Inclination < Values.Inclination + Tolerance.Inclination );

Bin.Velocity = ( Data.Velocity > Values.Velocity - Tolerance.Velocity ) & ...
               ( Data.Velocity < (Values.Velocity + Tolerance.Velocity) ); 

Bin.Slip.Angle = ( Data.Slip.Angle > Values.Slip.Angle - Tolerance.Slip.Angle ) & ...
                 ( Data.Slip.Angle < Values.Slip.Angle + Tolerance.Slip.Angle ); 

Bin.Slip.Ratio = ( Data.Slip.Ratio > Values.Slip.Ratio - Tolerance.Slip.Ratio ) & ...
                 ( Data.Slip.Ratio < Values.Slip.Ratio + Tolerance.Slip.Ratio ); 

% Eliminating Sparse Bins
for i = 1 : length( Values.Pressure )
    if sum( Bin.Pressure(i,:) ) < 50
        Bin.Pressure(i,:) = false( size( Bin.Pressure(i,:) ) );
    end   
end

for i = 1 : length( Values.Load )
    if sum( Bin.Load(i,:) ) < 50
        Bin.Load(i,:) = false( size( Bin.Load(i,:) ) );
    end    
end

for i = 1 : length( Values.Inclination )
    if sum( Bin.Inclination(i,:) ) < 50
        Bin.Inclination(i,:) = false( size( Bin.Inclination(i,:) ) );
    end    
end

for i = 1 : length( Values.Velocity )
    if sum( Bin.Velocity(i,:) ) < 50
         Bin.Velocity(i,:) = false( size( Bin.Velocity(i,:) ) );
    end    
end

Bin.Values = Values;
Bin.Tolerance = Tolerance;

end
