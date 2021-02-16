function Tire = PureLongitudinalFitting( Tire, Data, Bin, Figure )
%% PureLongitudinalFitting - Fits Pure Slip Longitudinal Force Pacejka Model 
% Executes all of the fitting procedures for pure slip longitudinal force
% generation. All equations are referenced from the 3rd Edition of 'Tyre &
% Vehicle Dynamics' by Pajecka.
%   The fitting process first fit the P6 form of the Pacejka Model and then 
%   fits a response surface to each of the primary parameters. The response
%   surface parameters are then used to initialize and bound the full
%   variant fit.
%
% Inputs:
%   Tire      - Tire Model
%   Data      - Parsed FSAE TTC Data
%   Bin       - Logical Binnings for Separating Operating Conditions
%   Figure    - Stores Model Figures
%
% Inputs:
%   Tire      - Tire Model w/ Pure Slip Longitudinal Force Model
%
% Author(s): 
% Blake Christierson (bechristierson@ucdavis.edu) [Sep 2018 - Jun 2021] 
% Carlos Lopez       (calopez@ucdavis.edu       ) [Jan 2019 -         ]
% 
% Last Updated: 15-Feb-2021

%% Operating Condition Space
Case.Pressure    = Bin(5).Values.Pressure;    % Pressure Bin Values Storage
Case.Load        = Bin(5).Values.Load;        % Normal Load Bin Values Storage
Case.Inclination = Bin(5).Values.Inclination; % Inclination Bin Values Storage

Mesh = struct( 'Pressure', [], 'Load', [], 'Inclination', [], 'dPi', [], 'dFz', [] );

for p = 1 : numel( Case.Pressure )
    for z = 1 : numel( Case.Load )
        for c = 1 : numel( Case.Inclination )
            Mesh(p,z,c).Pressure    = Case.Pressure(p);
            Mesh(p,z,c).Load        = Case.Load(z);
            Mesh(p,z,c).Inclination = Case.Inclination(c);
            
            Mesh(p,z,c).dPi = (Case.Pressure(p) - Tire.Pacejka.Pio) ./ Tire.Pacejka.Pio;
            Mesh(p,z,c).dFz = (Case.Load(z)     - Tire.Pacejka.Fzo) ./ Tire.Pacejka.Fzo;
        end
    end
end

%% Data Allocation
Raw = struct( 'Slip', [], 'Force'      , [], 'Pressure', [], ...
              'Load', [], 'Inclination', [], 'dFz', [], 'dPi', [] );
Raw( size(Mesh,1), size(Mesh,2), size(Mesh,3) ).Slip = [];

for i = [5 6]
    if isempty( Data(i).Source )
        continue
    end
    
    for p = 1 : numel( Case.Pressure )       
        for z = 1 : numel( Case.Load ) 
            for c = 1 : numel( Case.Inclination )
                Idx.Valid = Bin(i).Pressure(p,:) & Bin(i).Load(z,:) & Bin(i).Inclination(c,:) & ...
                    Bin(i).Slip.Angle( find( Bin(i).Values.Slip.Angle == 0 ), : );
                
                if sum( Idx.Valid ) < 50
                    continue % Skip Sparse Bins
                elseif (i == 6) && (Case.Pressure(p) == 12)
                    continue % Skip Tire Aging Sweep at 12 psi in Cornering 2
                end
                
                Raw(p,z,c).Slip  = Data(i).Slip.Ratio(Idx.Valid); % Allocate Slip Angle Data
                Raw(p,z,c).Force = Data(i).Force(1,Idx.Valid);    % Allocate Lateral Force Data
                
                Raw(p,z,c).Pressure    = Data(i).Pressure(Idx.Valid);    % Allocate Pressure Data
                Raw(p,z,c).Load        = Data(i).Force(3,Idx.Valid);     % Allocate Normal Force Data
                Raw(p,z,c).Inclination = Data(i).Inclination(Idx.Valid); % Allocate Inclination Data
                
                Raw(p,z,c).dFz = (Raw(p,z,c).Load     - Tire.Pacejka.Fzo) ./ Tire.Pacejka.Fzo;
                Raw(p,z,c).dPi = (Raw(p,z,c).Pressure - Tire.Pacejka.Pio) ./ Tire.Pacejka.Pio;
            end
        end
    end
end

%% Nominal Fitting
Nominal = struct( 'C0', NaN, 'D0', NaN, 'E0', NaN, 'K0', NaN, 'H0', NaN, 'V0', NaN, ...
                  'C' , NaN, 'D' , NaN, 'E' , NaN, 'K' , NaN, 'H' , NaN, 'V' , NaN, 'Residual', [] );
Nominal( size(Mesh,1), size(Mesh,2), size(Mesh,3) ).Residual = [];

for p = 1 : numel( Case.Pressure )
    for z = 1 : numel( Case.Load )
        for c = 1 : numel( Case.Inclination )
            if isempty( Raw(p,z,c).Slip )
                continue
            end
            
            Nominal(p,z,c) = PureLongitudinalNominal( Raw(p,z,c) );
        end
    end
end

%% Filtering Data & Operating Conditions
Mesh(    ind2sub(size(Raw), find(cellfun(@isempty, {Nominal.C}))) ) = [];
Raw(     ind2sub(size(Raw), find(cellfun(@isempty, {Nominal.C}))) ) = [];
Nominal(                         cellfun(@isempty, {Nominal.C})   ) = [];

Mesh(    ind2sub(size(Raw), find(cellfun(@isnan, {Nominal.C}))) ) = [];
Raw(     ind2sub(size(Raw), find(cellfun(@isnan, {Nominal.C}))) ) = [];
Nominal(                         cellfun(@isnan, {Nominal.C})   ) = [];

%% Variant Fitting
Response = PureLongitudinalResponseSurfaces( Tire, Mesh, Nominal );
  
[ Variant, Tire ] = PureLongitudinalVariant( Tire, Raw, Response );

%% Plotting Function
PureLongitudinalPlotting( Tire, Raw, Mesh, Nominal, Response, Variant, Figure );
