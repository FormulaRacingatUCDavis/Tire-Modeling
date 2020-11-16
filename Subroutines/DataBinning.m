function Bin = DataBinning( Data )

% Conversion Functions
lbf2N = @(lbf) lbf.*4.4482216152605;
mph2kmh = @(mph) mph.*1.6093440;

% Bin Values
Values.Pressure = [ 8 10 12 14 ]'; % Pressure Bins [psi]
Values.Load = lbf2N([ 50 100 150 200 250 350 ])'; % Normal Force Bins [N]
Values.Camber = [ -2 0 2 4 ]'; % Camber Bins [deg]
Values.Velocity = mph2kmh([ 0 2 15 25 45 ])'; % Belt Speed Bins [kph]
Values.Slip.Angle = [ -1 0 1 6]'; % Slip Angle Bins [deg]
Values.Slip.Ratio = 0; % Slip Ratio Bins [deg]

Tolerance.Pressure = min( diff(Values.Pressure) ) / 3; % Pressure Tolerance [psi]
Tolerance.Load = min( diff(Values.Load) ) / 3; % Normal Force Tolerance [N]
Tolerance.Camber = min( diff(Values.Camber) ) / 3; % Camber Tolerance [deg]
Tolerance.Velocity = min( diff(Values.Velocity) ) / 3; % Belt Velocity Tolerance [kph]
Tolerance.Slip.Angle = min( diff(Values.Slip.Angle) ) / 5; % Slip Angle Tolerance [deg]
Tolerance.Slip.Ratio = Tolerance.Slip.Angle * 0.15 / 12; % Slip Ratio Tolerance [deg]

% Bin Creation
Bin.Pressure = ( Data.Pressure > (Values.Pressure - Tolerance.Pressure) ) & ...
    ( Data.Pressure < (Values.Pressure + Tolerance.Pressure) );

Bin.Load = ( Data.Force(3,:) > (Values.Load - Tolerance.Load) ) & ...
    ( Data.Force(3,:) < (Values.Load + Tolerance.Load) );

Bin.Camber = ( Data.Camber > Values.Camber - Tolerance.Camber ) & ...
    ( Data.Camber < Values.Camber + Tolerance.Camber );

Bin.Velocity = ( Data.Velocity > Values.Velocity - Tolerance.Velocity ) & ...
    ( Data.Velocity < (Values.Velocity + Tolerance.Velocity) ); 

Bin.Slip.Angle = ( Data.Slip.Angle > Values.Slip.Angle - Tolerance.Slip.Angle ) & ...
    ( Data.Slip.Angle < Values.Slip.Angle + Tolerance.Slip.Angle ); 

Bin.Slip.Ratio = ( Data.Slip.Ratio > Values.Slip.Ratio - Tolerance.Slip.Ratio ) & ...
    ( Data.Slip.Ratio < Values.Slip.Ratio + Tolerance.Slip.Ratio ); 

Bin.Gain.Slip.Angle = ( diff( abs( Data.Slip.Angle ) ) >= 0);
Bin.Gain.Slip.Angle(end+1) = Bin.Gain.Slip.Angle(end); % Slip Angle Gain Bin

Bin.Gain.Slip.Ratio = ( diff( abs( Data.Slip.Ratio ) ) >= 0);
Bin.Gain.Slip.Ratio(end+1) = Bin.Gain.Slip.Ratio(end); % Slip Ratio Gain Bin

if strcmp( Data.TestName, 'Transient' )
    a = 1;
end

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

for i = 1 : length( Values.Camber )
    if sum( Bin.Camber(i,:) ) < 50
        Bin.Camber(i,:) = false( size( Bin.Camber(i,:) ) );
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
