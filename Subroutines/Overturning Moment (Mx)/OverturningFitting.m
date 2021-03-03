function Tire = OverturningFitting(Tire, Data, Bin, Figure)
%% Overturning Fitting = Allocates all the data for Variant and runs the plotting function
% Executes all the fitting procedures for the Overturning Moment
% generation. All equations are referenced from the 3rd Edition of 'Tyre &
% Vehicle Dynamics' (Page 182) from Pacejka.

%   The Fitting process fits the variant coefficient fitting at neutral
%   camber which is then used to fit camber variance. This will be modified
%   in the future for pressure variance.

% Fit Surface Variant Inclination & Pressure 
% Fit Primary Curves for Bounds & Initial Constrained High Dimensional 
% fmincon()
% Inputs:
%   Tire      - Tire Model
%   Data      - Parsed FSAE TTC Data
%   Bin       - Logical Binnings for Separating Operating Conditions
%   Figure    - Stores Model Figures
% Outputs:
%   Tire      - Tire Model w/ Relaxation Length Model
%
% Author(s): 
% Blake Christierson (bechristierson@ucdavis.edu) [Sep 2018 - Jun 2021]
% Carlos Lopez (calopez@ucdavis.edu) [Dec 2020 - ]
% 
% Last Updated: 20-Feb-2021
%% Operating Condition Space

Case.Pressure    = Bin(2).Values.Pressure;    % Pressure Bin Values Storage
Case.Load        = Bin(2).Values.Load;        % Normal Load Bin Values Storage
Case.Inclination = Bin(2).Values.Inclination; % Inclination Bin Values Storage

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
Raw = struct( 'Alpha', [], 'Kappa', [], 'Moment' , [], 'Pressure', [], ...
              'Load', [], 'Inclination', [], 'dFz', [], 'dPi', [] );
Raw( size(Mesh,1), size(Mesh,2), size(Mesh,3) ).Alpha = [];

for i = [2 3]
    if isempty( Data(i).Source )
        continue
    end

    for p = 1 : numel( Case.Pressure )       
        for z = 1 : numel( Case.Load ) 
            for c = 1 : numel( Case.Inclination )
                Idx.Valid = Bin(i).Pressure(p,:) & Bin(i).Load(z,:) & ...
                    Bin(i).Inclination(c,:); 
                
                if numel( Idx.Valid ) < 50
                    continue % Skip Sparse Bins
                elseif (i == 3) && (Case.Pressure(p) == 12)
                    continue % Skip Tire Aging Sweep at 12 psi in Cornering 2
                end
                
                Raw(p,z,c).Alpha  = Data(i).Slip.Angle(Idx.Valid); % Allocate Slip Angle Data
                Raw(p,z,c).Kappa  = Data(i).Slip.Ratio(Idx.Valid); % Allocate Slip Ratio Data
                Raw(p,z,c).Moment = Data(i).Moment(1, Idx.Valid);  % Allocate Overturning Moment Data
                
                Raw(p,z,c).Pressure    = Data(i).Pressure(Idx.Valid);    % Allocate Pressure Data
                Raw(p,z,c).Load        = Data(i).Force(3,Idx.Valid);     % Allocate Normal Force Data
                Raw(p,z,c).Inclination = Data(i).Inclination(Idx.Valid); % Allocate Inclination Data
                
                Raw(p,z,c).dFz = (Raw(p,z,c).Load     - Tire.Pacejka.Fzo) ./ Tire.Pacejka.Fzo;
                Raw(p,z,c).dPi = (Raw(p,z,c).Pressure - Tire.Pacejka.Pio) ./ Tire.Pacejka.Pio;
            end
        end
    end
end

%% Variant Fitting 
[Variant, Tire ] = OverturningVariant( Raw, Tire);

%% Plotting Function
OverturningPlotting( Mesh, Raw, Variant, Tire, Figure );

end