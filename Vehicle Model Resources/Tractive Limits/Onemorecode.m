% Define the vehicle parameters
vehicle_mass = 150;  % kg
wheel_radius = 0.3;  % m
g = 9.81;  % m/s^2
mu = 0.8;  % tire friction coefficient
Crr = 0.005;  % rolling resistance coefficient

% Define the velocity range
velocity_min = 0;  % m/s
velocity_max = 80;  % m/s
velocity_step = 10;  % m/s
velocity_range = velocity_min:velocity_step:velocity_max;

% Calculate the tractive limit at each velocity
tractive_limit = zeros(size(velocity_range));
for i = 1:length(velocity_range)
  velocity = velocity_range(i);
  tire_kappa = velocity / (wheel_radius * g); %Or can try to use function
  tractive_limit(i) = (mu * vehicle_mass * g * wheel_radius * (1 - tire_kappa)) / ...
    (1 + (Crr / mu)) - (vehicle_mass * g * Crr / mu);
end

% Plot the tractive limit vs. velocity
plot(velocity_range, tractive_limit);
xlabel('Velocity (m/s)');
ylabel('Tractive Limit (N)');