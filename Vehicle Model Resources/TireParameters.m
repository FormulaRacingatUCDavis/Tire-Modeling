function Tire = TireParameters(Name, Source, Notes)
%% TireParameters - Tire Model Structure Initialization
% This script generates an FRUCD tire model structure
% 
% Inputs:
%   Name   - (n,1 char)   Model Name
%   Source - (6,1 struct) Model Source File Information
%   Notes  - (n,1 cell)   Additional Model Notes
%
% Outputs:
%   Tire   - (1,1 struct) Initialized FRUCD Tire Model 
%
% Author(s): 
% Blake Christierson (bechristierson@ucdavis.edu) [Sep 2018 - Jun 2021] 
% 
% Last Updated: 02-May-2021

if nargin > 0
    Tire.Name    = Name;                      % Tire Name
    Tire.Date    = datestr(now,'dd-mm-yyyy'); % Generation Date
    Tire.Time    = datestr(now,'HH:MM PM'  ); % Generation Time
    Tire.Source  = Source;                    % Data Sources     
    Tire.Notes   = Notes;                     % Notes About Model

    Tire.Pacejka = PacejkaParameters();       % Pacejka Parameters
    Tire.Radius  = RadiusParameters();        % Radius Parameters
    Tire.Mass    = MassParameters();          % Mass Parameters
    Tire.Thermal = ThermalParameters();       % Thermal Parameters
end

function Pacejka = PacejkaParameters()
    %%% Nominal Conditions
    Pacejka.Fzo = 3000;       % Nominal Pacejka Load     [N] 
    Pacejka.Pio = 70;         % Nominal Pacejka Pressure [kPa]
    Pacejka.Ro  = 8*2.54/100; % Nominal Pacejka Radius   [m]
    Pacejka.Vo  = 15;         % Nominal Pacejka Velocity [m/s]
    
    % Scaling Factors (Lambdas)
    Pacejka.L.F.zo = 1; 

    Pacejka.L.mu.x = 1;
    Pacejka.L.mu.y = 1;
    Pacejka.L.mu.V = 0;

    Pacejka.L.K.x.k = 1;
    Pacejka.L.K.x.a = 1;

    Pacejka.L.K.y.a = 1;
    Pacejka.L.K.y.k = 1;

    Pacejka.L.K.y.g = 1;
    Pacejka.L.K.z.g = 1;

    Pacejka.L.C.x = 1;
    Pacejka.L.C.y = 1;
    Pacejka.L.C.z = 1;

    Pacejka.L.E.x = 1;
    Pacejka.L.E.y = 1;

    Pacejka.L.H.x = 1;
    Pacejka.L.H.y = 1;

    Pacejka.L.V.x = 1;
    Pacejka.L.V.y = 1;
    Pacejka.L.V.yk = 1;

    Pacejka.L.t = 1;

    Pacejka.L.M.x = 1;
    Pacejka.L.M.y = 1;
    Pacejka.L.M.r = 1;

    % Pacejka Spin Factors (Zetas)
    Pacejka.Z = ones(8,1);

    % p-Factors (Pure Slip Force Coefficients)
    Pacejka.p.C.x = 0;
    Pacejka.p.D.x = zeros(3,1);
    Pacejka.p.E.x = zeros(4,1);
    Pacejka.p.K.x = zeros(3,1);
    Pacejka.p.H.x = zeros(2,1);
    Pacejka.p.V.x = zeros(2,1);
    Pacejka.p.P.x = zeros(4,1);

    Pacejka.p.C.y = 0;
    Pacejka.p.D.y = zeros(3,1);
    Pacejka.p.E.y = zeros(5,1);
    Pacejka.p.K.y = zeros(7,1);
    Pacejka.p.H.y = zeros(2,1);
    Pacejka.p.V.y = zeros(4,1);
    Pacejka.p.P.y = zeros(5,1);

    % q-Factors (Moment Coefficients)
    Pacejka.q.B.z = zeros(10,1);
    Pacejka.q.C.z = 0;
    Pacejka.q.D.z = zeros(11,1);
    Pacejka.q.E.z = zeros(5,1);
    Pacejka.q.H.z = zeros(4,1);
    Pacejka.p.P.z = zeros(2,1);

    Pacejka.q.s.x = zeros(11,1);
    Pacejka.q.s.y = zeros(8,1);
    Pacejka.p.P.Mx = 0;
    
    % r-Factors (Combined Slip Force Coefficients)
    Pacejka.r.B.x = zeros(3,1);
    Pacejka.r.C.x = 0;
    Pacejka.r.E.x = zeros(2,1);
    Pacejka.r.H.x = 0;

    Pacejka.r.B.y = zeros(4,1);
    Pacejka.r.C.y = 0;
    Pacejka.r.E.y = zeros(2,1);
    Pacejka.r.H.y = zeros(2,1);
    Pacejka.r.V.y = zeros(6,1);

    % s-Factors (Combined Slip Aligning Moment Coefficients)
    Pacejka.s.s.z = zeros(4,1);
end

function Radius = RadiusParameters()
    Radius.Effective = [];
    Radius.Loaded    = [];
end

function Mass = MassParameters()
    Mass.m  = 0;            % Wheel Package Mass [kg]
    Mass.I  = zeros(3,3);   % Wheel Pacakge Inertia Tensor [kg-m^2]

    Mass.Is = 0;            % Spin Inertia (Spindle, Rotor, Rim, ...) [kg-m^2]
end

function Thermal = ThermalParameters()
    Thermal = [];
end

end

