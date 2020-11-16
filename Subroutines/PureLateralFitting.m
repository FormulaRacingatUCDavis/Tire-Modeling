function Tire = PureLateralFitting( Tire, Data, Bin )
% Executes all of the fitting procedures for pure slip lateral force
% generation. All equations are referenced from the 3rd Edition of 'Tyre &
% Vehicle Dynamics' by Pajecka.
%   The fitting process first fits the nominal coefficients to all cases. 
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

Mesh = struct( 'Pressure', [], 'Load', [], 'Camber', [], 'dPi', [], 'dFz', [] );

for p = 1 : numel( Case.Pressure )
    for z = 1 : numel( Case.Load )
        for c = 1 : numel( Case.Camber )
            Mesh(p,z,c).Pressure = Case.Pressure(p);
            Mesh(p,z,c).Load = Case.Load(z);
            Mesh(p,z,c).Camber = Case.Camber(c);
            
            Mesh(p,z,c).dPi = (Case.Pressure(p) - Tire.Pio) ./ Tire.Pio;
            Mesh(p,z,c).dFz = (Case.Load(z) - Tire.Fzo) ./ Tire.Fzo;
        end
    end
end

%% Data Allocation
Raw = struct( 'Slip', [], 'Force', [], 'Pressure', [], ...
    'Load', [], 'Camber', [], 'dFz', [], 'dPi', [] );
Raw( size(Mesh,1), size(Mesh,2), size(Mesh,3) ).Slip = [];

for i = [2 3]
    if isempty( Data(i).Source )
        continue
    end
    
    for p = 1 : numel( Case.Pressure )       
        for z = 1 : numel( Case.Load ) 
            for c = 1 : numel( Case.Camber )
                Idx.Valid = Bin(i).Pressure(p,:) & Bin(i).Load(z,:) & ...
                    Bin(i).Camber(c,:) & Bin(i).Gain.Slip.Angle & ...
                    Bin(i).Slip.Ratio( find( Bin(i).Values.Slip.Ratio == 0 ), : ); %#ok<FNDSB>
                
                if sum( Idx.Valid ) < 50
                    continue % SKip Sparse Bins
                elseif (i == 3) && (Case.Pressure(p) == 12)
                    continue % Skip Tire Aging Sweep at 12 psi in Cornering 2
                end
                
                Raw(p,z,c).Slip = Data(i).Slip.Angle(Idx.Valid); % Allocate Slip Angle Data
                Raw(p,z,c).Force = Data(i).Force(2,Idx.Valid); % Allocate Lateral Force Data
                
                Raw(p,z,c).Pressure = Data(i).Pressure(Idx.Valid); % Allocate Pressure Data
                Raw(p,z,c).Load = Data(i).Force(3,Idx.Valid); % Allocate Normal Force Data
                Raw(p,z,c).Camber = Data(i).Camber(Idx.Valid); % Allocate Camber Data
                
                Raw(p,z,c).dFz = (Raw(p,z,c).Load - Tire.Fzo) ./ Tire.Fzo;
                Raw(p,z,c).dPi = (Raw(p,z,c).Pressure - Tire.Pio) ./ Tire.Pio;
            end
        end
    end
end

%% Nominal Fitting
Nominal = struct( 'C0', NaN, 'D0', NaN, 'E0', NaN, 'K0', NaN, 'H0', NaN, 'V0', NaN, ...
    'C', NaN, 'D', NaN, 'E', NaN, 'K', NaN, 'H', NaN, 'V', NaN, 'Residual', [] );
Nominal( size(Mesh,1), size(Mesh,2), size(Mesh,3) ).Residual = [];

for p = 1 : numel( Case.Pressure )
    for z = 1 : numel( Case.Load )
        for c = 1 : numel( Case.Camber )
            if isempty( Raw(p,z,c).Slip )
                continue
            end
            
            Nominal(p,z,c) = PureLateralNominal( Raw(p,z,c) );
        end
    end
end

%% Filtering Data & Operating Conditions
for n = 1 : 3
    for i = size(Nominal,n): -1 : 1
        switch n
            case 1
                if isempty( [Nominal(i,:,:).Residual] )
                    Case.Pressure(i) = [];
                    Mesh(i,:,:) = [];
                    Raw(i,:,:) = [];
                    Nominal(i,:,:) = [];
                end
            case 2
                if isempty( [Nominal(:,i,:).Residual] )
                    Case.Load(i) = [];
                    Mesh(:,i,:) = [];
                    Raw(:,i,:) = [];
                    Nominal(:,i,:) = [];
                end
            case 3
                if isempty( [Nominal(:,:,i).Residual] )
                    Case.Camber(i) = [];
                    Mesh(:,:,i) = [];
                    Raw(:,:,i) = [];
                    Nominal(:,:,i) = [];
                end
        end
    end
end

%% Variant Fitting
Response = PureLateralResponseSurfaces( Mesh, Nominal, Tire );
  
[ Variant, Tire ] = PureLateralVariant( Raw, Response.x0, Tire ); %#ok<ASGLU>

%% Plotting Function
PureLateralPlotting( Mesh, Raw, Nominal, Tire );
