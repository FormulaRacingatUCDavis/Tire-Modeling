% Define vehicle parameters
vehicleMass = 454; % kg
vehicleCGHeight = 0.55; % m
vehicleTrackWidth = 1.22; % m
vehicleWheelbase = 1.6; % m
vehicleFrontalArea = 0.7; % m^2
vehicleDragCoefficient = 0.3;

% Define tire parameters
tireRadius = 0.29; % m
tirePeakFrictionCoefficient = 1.2;

% Define acceleration and deceleration parameters
acceleration = 5; % m/s^2
braking = -5; % m/s^2

% Define speed range
speedRange = 0:10:100; % km/h

% Convert speed range to m/s
speedRange = speedRange / 3.6; % m/s

% Calculate aerodynamic drag force
dragForce = vehicleDragCoefficient * vehicleFrontalArea * speedRange.^2; % N

% Calculate rolling resistance force
rollingResistanceForce = 0.01 * vehicleMass * 9.81; % N

% Calculate weight transfer due to acceleration
weightTransfer = vehicleMass * acceleration * vehicleCGHeight / vehicleWheelbase; % N

% Calculate weight transfer due to braking
weightTransfer = vehicleMass * braking * vehicleCGHeight / vehicleWheelbase; % N

% Calculate lateral tire force
lateralTireForce = tirePeakFrictionCoefficient * vehicleMass * 9.81 / tireRadius; % N

% Calculate tractive limit
tractiveLimit = (dragForce + rollingResistanceForce + weightTransfer) / lateralTireForce; % dimensionless

% Plot tractive limit vs speed
plot(speedRange, tractiveLimit);
xlabel('Speed (m/s)');
ylabel('Tractive limit');