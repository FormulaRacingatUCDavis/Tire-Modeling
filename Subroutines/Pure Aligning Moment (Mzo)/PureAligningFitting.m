function Tire = PureAligningFitting( Tire, Data, Bin, Figure )
% Executes all of the fitting procedures for pure slip aligning moment
% generation. All equations are referenced from the 3rd Edition of 'Tyre &
% Vehicle Dynamics' by Pajecka.
%   The fitting process first fits the nominal coefficients to all cases. 
%   Statistical analysis is then done to set upper and lower bounds of 
%   these coefficients during full load variance fitting at neutral camber. 
%   These are then used to fit camber variance. In the future this may be
%   modified for pressure v ariance. 

% Nominal Fit for Primaries ( Bt, Ct, Dt, Et, Ht, Br, Cr, Dr, Hf )
% Fit Surface Variant Inclination & Pressure
% Fit Primary Curves for Bounds & Initial
% Constrained High Dimensional fmincon()

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
Raw = struct( 'Slip'    , [], 'Force', [], 'Moment'     , [], ...
              'Pressure', [], 'Load' , [], 'Inclination', [], 'dFz', [], 'dPi', [] );
Raw( size(Mesh,1), size(Mesh,2), size(Mesh,3) ).Slip = [];

for i = [2 3]
    if isempty( Data(i).Source )
        continue
    end
    
    for p = 1 : numel( Case.Pressure )       
        for z = 1 : numel( Case.Load ) 
            for c = 1 : numel( Case.Inclination )
                Idx.Valid = Bin(i).Pressure(p,:) & Bin(i).Load(z,:) & Bin(i).Inclination(c,:) & ...
                    Bin(i).Slip.Ratio( find( Bin(i).Values.Slip.Ratio == 0 ), : ) & ...
                    ( abs(Data(i).Slip.Angle) < deg2rad(11) ); 
                
                if sum( Idx.Valid ) < 50
                    continue % Skip Sparse Bins
                elseif (i == 3) && (Case.Pressure(p) == 12)
                    continue % Skip Tire Aging Sweep at 12 psi in Cornering 2
                end
                
                Raw(p,z,c).Slip   = Data(i).Slip.Angle(Idx.Valid); % Allocate Slip Angle Data
                Raw(p,z,c).Force  = Data(i).Force(2,Idx.Valid);    % Allocate Lateral Force Data
                Raw(p,z,c).Moment = Data(i).Moment(3,Idx.Valid);   % Allocate Aligning Moment Data 
                
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
Nominal = struct( 'Bt'   , NaN, 'Ct', NaN, 'Dt', NaN, 'Et', NaN, 'Ht', NaN, ...
                  'qbz10', NaN, 'Cr', 1  , 'Dr', NaN, 'Hf', NaN, 'Residual', [] );
Nominal( size(Mesh,1), size(Mesh,2), size(Mesh,3) ).Residual = [];

for p = 1 : numel( Case.Pressure )
    for z = 1 : numel( Case.Load )
        for c = 1 : numel( Case.Inclination )
            if isempty( Raw(p,z,c).Slip )
                continue
            end
            
            Nominal(p,z,c) = PureAligningNominal( Tire, Raw(p,z,c), Mesh(p,z,c) );
        end
    end
end

%% Filtering Data & Operating Conditions
Mesh(    ind2sub(size(Raw), find(cellfun(@isempty, {Nominal.Ct}))) ) = [];
Raw(     ind2sub(size(Raw), find(cellfun(@isempty, {Nominal.Ct}))) ) = [];
Nominal(                         cellfun(@isempty, {Nominal.Ct})   ) = [];

Mesh(    ind2sub(size(Raw), find(cellfun(@isnan, {Nominal.Ct}))) ) = [];
Raw(     ind2sub(size(Raw), find(cellfun(@isnan, {Nominal.Ct}))) ) = [];
Nominal(                         cellfun(@isnan, {Nominal.Ct})   ) = [];

%% Variant Fitting
Response = PureAligningResponseSurfaces( Tire, Mesh, Nominal );
  
[ Variant, Tire ] = PureAligningVariant( Tire, Raw, Nominal, Response );

%% Plotting Function
PureAligningPlotting( Tire, Raw, Mesh, Nominal, Response, Variant, Figure );
