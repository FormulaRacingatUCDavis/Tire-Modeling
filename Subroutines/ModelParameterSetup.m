function Tire = ModelParameterSetup( )
    %% Pacejka Coefficients
    %%% Nominal Conditions
    Tire.Fzo = 2000    ; % Nominal Tire Load [N] (Chosen as 450lbs)
    Tire.Pio = 10      ; % Nominal Tire Pressure [psi]
    Tire.Ro  = 8 * 25.4; % Nominal Tire Radius [m]

    %%% Scaling Factors (Lambdas)
    Tire.L.F.zo = 1; 

    Tire.L.mu.x = 1;
    Tire.L.mu.y = 1;
    Tire.L.mu.V = 0;

    Tire.L.K.x.k = 1;
    Tire.L.K.x.a = 1;

    Tire.L.K.y.a = 1;
    Tire.L.K.y.k = 1;

    Tire.L.K.y.g = 1;
    Tire.L.K.z.g = 1;

    Tire.L.C.x = 1;
    Tire.L.C.y = 1;
    Tire.L.C.z = 1;

    Tire.L.E.x = 1;
    Tire.L.E.y = 1;

    Tire.L.H.x = 1;
    Tire.L.H.y = 1;

    Tire.L.V.x = 1;
    Tire.L.V.y = 1;
    Tire.L.V.yk = 1;

    Tire.L.t = 1;

    Tire.L.M.x = 1;
    Tire.L.M.y = 1;
    Tire.L.M.r = 1;

    %%% Tire Spin Factors (Zetas)
    Tire.Z = ones(8,1);

    %%% p-Factors
    Tire.p.C.x = 0;
    Tire.p.D.x = zeros(3,1);
    Tire.p.E.x = zeros(4,1);
    Tire.p.K.x = zeros(3,1);
    Tire.p.H.x = zeros(2,1);
    Tire.p.V.x = zeros(2,1);
    Tire.p.P.x = zeros(4,1);

    Tire.p.C.y = 0;
    Tire.p.D.y = zeros(3,1);
    Tire.p.E.y = zeros(5,1);
    Tire.p.K.y = zeros(7,1);
    Tire.p.H.y = zeros(2,1);
    Tire.p.V.y = zeros(4,1);
    Tire.p.P.y = zeros(5,1);

    %%% q-Factors
    Tire.q.B.z = zeros(10,1);
    Tire.q.C.z = 0;
    Tire.q.D.z = zeros(11,1);
    Tire.q.E.z = zeros(5,1);
    Tire.q.H.z = zeros(4,1);

    %%% r-Factors (Combined Slip)
    Tire.r.B.x = zeros(3,1);
    Tire.r.C.x = 0;
    Tire.r.H.x = 0;

    Tire.r.B.y = zeros(4,1);
    Tire.r.C.y = 0;
    Tire.r.H.y = 0;
    Tire.r.V.y = zeros(6,1);

    %%% s-Factors
    Tire.s.s.z = zeros(4,1);

    %% Radial Factors
    Tire.Re = [];
    Tire.Rl = [];
end