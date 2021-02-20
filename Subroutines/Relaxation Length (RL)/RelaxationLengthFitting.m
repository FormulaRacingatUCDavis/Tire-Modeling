function Tire = RelaxationLengthFitting( Tire, Data, Bin, Figure )
%% RelaxationLengthFitting - Fits Relaxation Length Pacejka Model 
% Executes all of the fitting procedures for relaxation length
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
% Outputs:
%   Tire      - Tire Model w/ Relaxation Length Model
%
% Author(s): 
% Blake Christierson (bechristierson@ucdavis.edu) [Sep 2018 - Jun 2021] 
% Leonardo Howard    (leohoward@ucdavis.edu     ) [Feb 2021 -         ]
% 
% Last Updated: 17-Feb-2021

%% Operating Condition Space
Case.Pressure    = Bin(1).Values.Pressure;    % Pressure Bin Values Storage
Case.Load        = Bin(1).Values.Load;        % Normal Load Bin Values Storage

Mesh = struct( 'Pressure', [], 'Load', [], 'dPi', [], 'dFz', [] );

for p = 1 : numel( Case.Pressure )
    for z = 1 : numel( Case.Load )
        Mesh(p,z).Pressure    = Case.Pressure(p);
        Mesh(p,z).Load        = Case.Load(z);
                        
        Mesh(p,z).dPi = (Case.Pressure(p) - Tire.Pacejka.Pio) ./ Tire.Pacejka.Pio;
        Mesh(p,z).dFz = (Case.Load(z)     - Tire.Pacejka.Fzo) ./ Tire.Pacejka.Fzo; 
    end
end

%% Data Allocation
Raw = struct( 'Slip', [], 'Force', [], 'Velocity', [], 'Pressure', [], ...
              'Load', [], 'dFz', [], 'dPi', [] );
Raw( size(Mesh,1), size(Mesh,2) ).Slip = [];

for p = 1 : numel( Case.Pressure )       
    for z = 1 : numel( Case.Load ) 
        Idx.Valid = find(Bin(1).Pressure(p,:) & Bin(1).Load(z,:));
        Idx.Valid = Idx.Valid( Idx.Valid > find( Bin(1).Velocity(4,:), 1, 'last' ) );
        % finds last point of initial spring rate tests

        if numel( Idx.Valid ) < 50
            continue % Skip Sparse Bins
        end

        Raw(p,z).Slip = Data(1).Slip.Angle(Idx.Valid);    % Allocate Slip Angle Data
        Raw(p,z).Force = Data(1).Force(2,Idx.Valid);    % Allocate Lateral Force Data
        Raw(p,z).Velocity = Data(1).Velocity(Idx.Valid);    % Allocate Velocity Data

        Raw(p,z).Pressure    = Data(1).Pressure(Idx.Valid);    % Allocate Pressure Data
        Raw(p,z).Load        = Data(1).Force(3,Idx.Valid);     % Allocate Normal Force Data

        Raw(p,z).dFz = (Raw(p,z).Load     - Tire.Pacejka.Fzo) ./ Tire.Pacejka.Fzo;
        Raw(p,z).dPi = (Raw(p,z).Pressure - Tire.Pacejka.Pio) ./ Tire.Pacejka.Pio;
    end
end

%% Filtering Data & Operating Conditions
Mesh(    ind2sub(size(Raw), find(cellfun(@isempty, {Raw.dFz}))) ) = [];
Raw(     ind2sub(size(Raw), find(cellfun(@isempty, {Raw.dFz}))) ) = [];

Mesh(    ind2sub(size(Raw), find(cellfun(@isnan, {Raw.dFz}))) ) = [];
Raw(     ind2sub(size(Raw), find(cellfun(@isnan, {Raw.dFz}))) ) = [];

%% Variant Fitting
Response = RelaxationLengthResponseSurfaces( Tire, Mesh, Nominal );
  
[ Variant, Tire ] = RelaxationLengthVariant( Tire, Raw, Response );

%% Plotting Function
RelaxationLengthPlotting( Tire, Raw, Mesh, Nominal, Response, Variant, Figure );
