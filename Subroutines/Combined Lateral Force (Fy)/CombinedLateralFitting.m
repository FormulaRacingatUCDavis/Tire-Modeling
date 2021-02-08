function Tire = CombinedLateralFitting(Tire, Data, Bin)
%Executes all fittings procedures for combined slip lateral force 
%generation. All equations are referenced from the 3rd Edition 
%(Section 4.3.2) of "Tyre and Vehicle Dynamics" by Pacejka
%   Fitting process first fits nominal coefficients to all cases.
%   Statistical analysis is then done to set upper and lower bounds of 
%   these coefficients during full load variance fitting at neutral camber. 
%   These are then used to fit camber variance. In the future this may be
%   modified for pressure variance. 

% Nominal Fit for Primaries (C,D,E,K,S_H,S_V)
% Fit Surface Variant Camber & Pressure
% Fit Primary Curves for Bounds & Initial
% Constrained High Dimensional fmincon()
%% Operating Condition Space
Case.Pressure = Bin(2).Values.Pressure; % Pressure Bin Values Storage
Case.Load = Bin(2).Values.Load; % Normal Load Bin Values Storage
Case.Camber = Bin(2).Values.Camber; % Camber Bin Values Storage
Case.Slip = Bin(2).Values.Slip; %Slip Bin Values Storage

Mesh = struct( "Pressure", [], "Load", [], "Camber", [],"SA0", [], ...
     "SA3", [], "SA6", [], "dPi", [],"dFz", [] );
%% Data Allocation
Raw = struct( "Slip", [], "Force", [], "Pressure", [], ...
    "Load", [], "Camber", [], "dFz", [], "dPi", [] );
%% Nominal Fitting
Nominal = struct( 'C0', NaN, 'D0', NaN, 'E0', NaN, 'K0', NaN, 'H0', NaN, 'V0', NaN, ...
    'C', NaN, 'D', NaN, 'E', NaN, 'K', NaN, 'H', NaN, 'V', NaN, 'Residual', [] );
Nominal( size(Mesh,1), size(Mesh,2), size(Mesh,3), size(Mesh,4)).Residual = [];
%% Filthering Data & Operating Conditions
for n = 1 : 4
    for i = size(Nominal,n): -1 : 1
        switch n
            case 1
                if isempty( [Nominal(i,:,:,:).Residual] )
                    Case.Pressure(i) = [];
                    Mesh(i,:,:) = [];
                    Raw(i,:,:) = [];
                    Nominal(i,:,:) = [];
                end
            case 2
                if isempty( [Nominal(:,i,:,:).Residual] )
                    Case.Load(i) = [];
                    Mesh(:,i,:) = [];
                    Raw(:,i,:) = [];
                    Nominal(:,i,:) = [];
                end
            case 3
                if isempty( [Nominal(:,:,i,:).Residual] )
                    Case.Camber(i) = [];
                    Mesh(:,:,i) = [];
                    Raw(:,:,i) = [];
                    Nominal(:,:,i) = [];
                end
            case 4
                if isempty( [Nominal(:,:,:,i).Residual])
                    Case.SA(i) = [];
                    Mesh(:,:,:,i) = [];
                    Raw(:,:,:,i) = [];
                    Nominal(:,:,:,i) = [];
                end
                    
        end
    end
end

%% Variant Fitting
Response = CombinedLateralResponseSurfaces(Mesh, Nominal, Tire);

[Variant, Tire] = CombinedLateralVariant(Raw, Response.x0, Tire);
%% Plotting Function
CombinedLateralPlotting(Mesh, Raw, Nominal, Tire);
