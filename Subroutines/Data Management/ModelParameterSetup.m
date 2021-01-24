function Tire = ModelParameterSetup( )
% Allocates model Coefficients into Structure

%% Pacejka Coefficients
% Nominal Conditions
Tire.Pacejka.Fzo = 3000; % Nominal Tire.Pacejka Load [N] (Chosen as 450lbs)
Tire.Pacejka.Pio = 10; % Nominal Tire.Pacejka Pressure [psi]
Tire.Pacejka.Ro = 8*2.54/100; % Nominal Tire.Pacejka Radius [m]

% Scaling Factors (Lambdas)
Tire.Pacejka.L.F.zo = 1; 

Tire.Pacejka.L.mu.x = 1;
Tire.Pacejka.L.mu.y = 1;
Tire.Pacejka.L.mu.V = 0;

Tire.Pacejka.L.K.x.k = 1;
Tire.Pacejka.L.K.x.a = 1;

Tire.Pacejka.L.K.y.a = 1;
Tire.Pacejka.L.K.y.k = 1;

Tire.Pacejka.L.K.y.g = 1;
Tire.Pacejka.L.K.z.g = 1;

Tire.Pacejka.L.C.x = 1;
Tire.Pacejka.L.C.y = 1;
Tire.Pacejka.L.C.z = 1;

Tire.Pacejka.L.E.x = 1;
Tire.Pacejka.L.E.y = 1;

Tire.Pacejka.L.H.x = 1;
Tire.Pacejka.L.H.y = 1;

Tire.Pacejka.L.V.x = 1;
Tire.Pacejka.L.V.y = 1;
Tire.Pacejka.L.V.yk = 1;

Tire.Pacejka.L.t = 1;

Tire.Pacejka.L.M.x = 1;
Tire.Pacejka.L.M.y = 1;
Tire.Pacejka.L.M.r = 1;

% Tire.Pacejka Spin Factors (Zetas)
Tire.Pacejka.Z = ones(8,1);

% p-Factors (Pure Slip Force Coefficients)
Tire.Pacejka.p.C.x = 0;
Tire.Pacejka.p.D.x = zeros(3,1);
Tire.Pacejka.p.E.x = zeros(4,1);
Tire.Pacejka.p.K.x = zeros(3,1);
Tire.Pacejka.p.H.x = zeros(2,1);
Tire.Pacejka.p.V.x = zeros(2,1);
Tire.Pacejka.p.P.x = zeros(4,1);

Tire.Pacejka.p.C.y = 0;
Tire.Pacejka.p.D.y = zeros(3,1);
Tire.Pacejka.p.E.y = zeros(5,1);
Tire.Pacejka.p.K.y = zeros(7,1);
Tire.Pacejka.p.H.y = zeros(2,1);
Tire.Pacejka.p.V.y = zeros(4,1);
Tire.Pacejka.p.P.y = zeros(5,1);

% q-Factors (Moment Coefficients)
Tire.Pacejka.q.B.z = zeros(10,1);
Tire.Pacejka.q.C.z = 0;
Tire.Pacejka.q.D.z = zeros(11,1);
Tire.Pacejka.q.E.z = zeros(5,1);
Tire.Pacejka.q.H.z = zeros(4,1);
Tire.Pacejka.p.P.z = zeros(2,1);

Tire.Pacejka.q.s.x = zeros(10,1);
Tire.Pacejka.q.s.y = zeros(8,1);

% r-Factors (Combined Slip Force Coefficients)
Tire.Pacejka.r.B.x = zeros(3,1);
Tire.Pacejka.r.C.x = 0;
Tire.Pacejka.r.E.x = zeros(2,1);
Tire.Pacejka.r.H.x = 0;

Tire.Pacejka.r.B.y = zeros(4,1);
Tire.Pacejka.r.C.y = 0;
Tire.Pacejka.r.E.y = zeros(2,1);
Tire.Pacejka.r.H.y = zeros(2,1);
Tire.Pacejka.Pacejka.r.V.y = zeros(6,1);

% s-Factors (Combined Slip Aligning Moment Coefficients)
Tire.Pacejka.s.s.z = zeros(4,1);

%% Radial Deflection Models
Tire.Radius.Effective = [];
Tire.Radius.Loaded    = [];

end