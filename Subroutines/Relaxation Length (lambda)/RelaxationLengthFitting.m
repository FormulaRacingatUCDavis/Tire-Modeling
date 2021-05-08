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
% Last Updated: 16-Apr-2021

%% Operating Condition Space
Case.Pressure    = Bin(1).Values.Pressure;    % Pressure Bin Values Storage
Case.Load        = Bin(1).Values.Load;        % Normal Load Bin Values Storage
Case.Slip.Angle  = Bin(1).Values.Slip.Angle;  % Slip Angle Bin Values Storage

Mesh = struct( 'Pressure', [], 'Load', [], 'dPi', [], 'dFz', [] );

%for a = 1 : numel( Case.Slip.Angle )
    for p = 1 : numel( Case.Pressure )
        for z = 1 : numel( Case.Load )
            Mesh(p,z).Pressure    = Case.Pressure(p);
            Mesh(p,z).Load        = Case.Load(z);
            
            Mesh(p,z).Slip.Angle  = Case.Slip.Angle;
            %Mesh(p,z).Slip.Angle  = Case.Slip.Angle(a);
            Mesh(p,z).Slip.Angle(2) = [];

            Mesh(p,z).dPi = (Case.Pressure(p) - Tire.Pacejka.Pio) ./ Tire.Pacejka.Pio;
            Mesh(p,z).dFz = (Case.Load(z)     - Tire.Pacejka.Fzo) ./ Tire.Pacejka.Fzo; 
        end
    end
%end

%% Data Allocation
Raw = struct( 'Slip', [], 'Force', [], 'Velocity', [], 'Pressure', [], ...
              'Load', [], 'dFz', [], 'dPi', [] );
Raw( size(Mesh,1), size(Mesh,2) ).Slip = [];

for p = 1 : numel( Case.Pressure )       
    for z = 1 : numel( Case.Load ) 
        Idx.Valid = find(Bin(1).Pressure(p,:) & Bin(1).Load(z,:));
        Idx.Valid = Idx.Valid( Idx.Valid > find( Bin(1).Velocity(4,:), 1, 'last' ) );     % Finds last point of initial spring rate tests
        
        if numel( Idx.Valid ) < 50
            continue % Skip Sparse Bins
        end
        
        Raw(p,z).Time       = Data(1).Time(Idx.Valid);          % Allocate Time Data
        Raw(p,z).Slip       = Data(1).Slip.Angle(Idx.Valid);    % Allocate Slip Angle Data
        Raw(p,z).Force      = Data(1).Force(2,Idx.Valid);       % Allocate Lateral Force Data
        Raw(p,z).Moment     = Data(1).Moment(3,Idx.Valid);      % Allocate Aligning Moment Data
        
        Raw(p,z).Velocity   = Data(1).Velocity(Idx.Valid);      % Allocate Velocity Data
        Raw(p,z).Load       = Data(1).Force(3,Idx.Valid);       % Allocate Normal Force Data
        Raw(p,z).Pressure   = Data(1).Pressure(Idx.Valid);      % Allocate Pressure Data
        
        Raw(p,z).dFz        = (Raw(p,z).Load     - Tire.Pacejka.Fzo) ./ Tire.Pacejka.Fzo;
        Raw(p,z).dPi        = (Raw(p,z).Pressure - Tire.Pacejka.Pio) ./ Tire.Pacejka.Pio;
    end
end

%% Filtering Data & Operating Conditions
Mesh(    ind2sub(size(Raw), find(cellfun(@isempty, {Raw.dFz}))) ) = [];
Raw(     ind2sub(size(Raw), find(cellfun(@isempty, {Raw.dFz}))) ) = [];

%% Velocity Thresholding
for j = 1 : numel(Raw)
    Idx.IsRolling = find( Raw(j).Velocity > ( Bin(1).Tolerance.Velocity ./ 3 ) ); 
    
%% Segmenting Velocity Based on Peak Locations
    [~, Idx.Breaks] = ...
        findpeaks( diff(Idx.IsRolling), 'NPeaks', 6, ...
        'MinPeakDistance', ( length( Raw(j).Velocity(Idx.IsRolling) ) ./ 10 ) );
   
    Idx.NewBreaks = Idx.Breaks( diff(Idx.Breaks) < mean( diff(Idx.Breaks) ) );
    
    for k = 1 : numel(Idx.NewBreaks)  
        Response(j,k).Idx = Idx.IsRolling( Idx.NewBreaks(k)+1 : ...
            Idx.Breaks( find( Idx.Breaks == Idx.NewBreaks(k) )+1 ) );
        
        Response(j,k).Time      = Raw(j).Time( Response(j,k).Idx );      % Allocate Transient Time Data
        Response(j,k).Slip      = Raw(j).Slip( Response(j,k).Idx );      % Allocate Transient Slip Angle Data
        Response(j,k).Force     = Raw(j).Force( Response(j,k).Idx );     % Allocate Transient Force Data
        Response(j,k).Moment    = Raw(j).Moment( Response(j,k).Idx );    % Allocate Transient Moment Data
        
        Response(j,k).Velocity  = Raw(j).Velocity( Response(j,k).Idx );  % Allocate Transient Velocity Data
        Response(j,k).Load      = Raw(j).Load( Response(j,k).Idx );      % Allocate Transient Normal Force Data
        Response(j,k).Pressure  = Raw(j).Pressure( Response(j,k).Idx );  % Allocate Transient Pressure Data
    end
end

%% Plotting Function
RelaxationLengthPlotting( Tire, Raw, Response, Idx, Mesh, Figure );