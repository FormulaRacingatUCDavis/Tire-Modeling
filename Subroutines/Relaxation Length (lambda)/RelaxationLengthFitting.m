function Tire = RelaxationLengthFitting( Tire, Data, Bin, Figure )
%% Relaxation Length Modeling - Transient Tire Model
% This fits surfaces for the tire relaxation length as a function of 
% inflation pressure and normal load for different slip angles. Raw data
% is segmented to where transient states occur to obtain variable tau from
% the curve fitting, which is used to calculate for the relaxation length.
%  
% Inputs:
%   Tire           - (struct)      Pacejka Parameters
%   SlipAngle      - (n,1 numeric) Slip Angle          {alpha} [deg]
%   NormalLoad     - (n,1 numeric) Normal Load         {F_z}   [N]
%   LateralForce   - (n,2 numeric) Lateral Force       {F_y}   [N]
%   AligningMoment - (n,3 numeric) Aligning Moment     {M_z}   [Nm]
%   Pressure       - (n,1 numeric) Inflation Pressure  {P_i}   [kPa]
%   Velocity       - (n,1 numeric) Center Velocity     {v_c}   [m/s]
%   Model          - (struct)      Fidelity Choices 
%
% Outputs:
%   RelaxationLength - (n,1 numeric) Relaxation Length {lambda} [m]
%
% Notes:
%   2nd order curve fitting is currently not implemented (6/2/21) because
%   accuracy of relaxation length data is currently sufficient to justify 
%   tire choice.
%
% Author(s): 
% Blake Christierson (bechristierson@ucdavis.edu) [Sep 2018 - Jun 2021] 
% Leonardo Howard    (leohoward@ucdavis.edu     ) [Feb 2021 -         ]
% 
% Last Updated: 02-June-2021

%% Operating Condition Space
for x = 1 : numel( Data )
    Case.Pressure    = Bin(x).Values.Pressure;    % Pressure Bin Values Storage
    Case.Load        = Bin(x).Values.Load;        % Normal Load Bin Values Storage
    Case.Slip.Angle  = Bin(x).Values.Slip.Angle;  % Slip Angle Bin Values Storage

    Mesh = struct( 'Pressure', [], 'Load', [], 'dPi', [], 'dFz', [] );

    for p = 1 : numel( Case.Pressure )
        for z = 1 : numel( Case.Load )
            Mesh(p,z).Pressure    = Case.Pressure(p);
            Mesh(p,z).Load        = Case.Load(z);

            Mesh(p,z).Slip.Angle  = Case.Slip.Angle;
            Mesh(p,z).Slip.Angle(2) = [];

            Mesh(p,z).dPi = (Case.Pressure(p) - Tire(x).Pacejka.Pio) ./ Tire(x).Pacejka.Pio;
            Mesh(p,z).dFz = (Case.Load(z)     - Tire(x).Pacejka.Fzo) ./ Tire(x).Pacejka.Fzo; 
        end
    end

%% Data Allocation
    Raw = struct( 'Slip', [], 'Force', [], 'Velocity', [], 'Pressure', [], ...
                  'Load', [], 'dFz', [], 'dPi', [] );
    Raw( size(Mesh,1), size(Mesh,2) ).Slip = [];

    for p = 1 : numel( Case.Pressure )       
        for z = 1 : numel( Case.Load ) 
            Idx.Valid = find(Bin(x).Pressure(p,:) & Bin(x).Load(z,:));

            % Finds last point of initial spring rate tests
            Idx.Valid = Idx.Valid( Idx.Valid > find( Bin(x).Velocity(4,:), 1, 'last' ) );     

            if numel( Idx.Valid ) < 50
                continue % Skip Sparse Bins
            end

            Raw(p,z).Time       = Data(x).Time(Idx.Valid);          % Allocate Time Data
            Raw(p,z).Slip       = Data(x).Slip.Angle(Idx.Valid);    % Allocate Slip Angle Data
            Raw(p,z).Force      = Data(x).Force(2,Idx.Valid);       % Allocate Lateral Force Data
            Raw(p,z).Moment     = Data(x).Moment(3,Idx.Valid);      % Allocate Aligning Moment Data

            Raw(p,z).Velocity   = Data(x).Velocity(Idx.Valid);      % Allocate Velocity Data
            Raw(p,z).Load       = Data(x).Force(3,Idx.Valid);       % Allocate Normal Force Data
            Raw(p,z).Pressure   = Data(x).Pressure(Idx.Valid);      % Allocate Pressure Data

            Raw(p,z).dFz        = (Raw(p,z).Load     - Tire(x).Pacejka.Fzo) ./ Tire(x).Pacejka.Fzo;
            Raw(p,z).dPi        = (Raw(p,z).Pressure - Tire(x).Pacejka.Pio) ./ Tire(x).Pacejka.Pio;
        end
    end
%% Filtering Data & Operating Conditions
    Mesh(    ind2sub(size(Raw), find(cellfun(@isempty, {Raw.dFz}))) ) = [];
    Raw(     ind2sub(size(Raw), find(cellfun(@isempty, {Raw.dFz}))) ) = [];

%% Rising Edge Detection
    %{ 
        To segment out the data where transient states occur, raw tire 
        data is filtered based on when tire is rolling. Then, the data is 
        further filtered based on velocity peaks because transient states 
        occur when the velocity reaches maximum and there is a total of 6 
        instances per case, which is eventually trimmed down to 3 instances
        per case because relaxation length is in regard to the delay in
        the buildup of the cornering force and it occurs only when the
        number of data points in between peaks is less than its average.
    %}

    PressureCase=1;
    SlipCase=1;
    
    for j = 1 : numel(Raw)
        Idx.IsRolling = find( Raw(j).Velocity > ( Bin(1).Tolerance.Velocity ./ 3 ) ); 

        [~, Idx.Breaks] = ...
            findpeaks( diff(Idx.IsRolling), 'NPeaks', 6, ...
            'MinPeakDistance', ( length( Raw(j).Velocity(Idx.IsRolling) ) ./ 10 ) );

        Idx.NewBreaks = Idx.Breaks( diff(Idx.Breaks) < mean( diff(Idx.Breaks) ) );

        for k = 1 : numel(Idx.NewBreaks)  
            Run(x).Response(j,k).Idx = Idx.IsRolling( Idx.NewBreaks(k)+1 : ...
                Idx.Breaks( find( Idx.Breaks == Idx.NewBreaks(k) )+1 ) );

            Run(x).Response(j,k).Time      = Raw(j).Time( Run(x).Response(j,k).Idx );      % Allocate Transient Time Data
            Run(x).Response(j,k).Slip      = Raw(j).Slip( Run(x).Response(j,k).Idx );      % Allocate Transient Slip Angle Data
            Run(x).Response(j,k).Force     = Raw(j).Force( Run(x).Response(j,k).Idx );     % Allocate Transient Force Data
            Run(x).Response(j,k).Moment    = Raw(j).Moment( Run(x).Response(j,k).Idx );    % Allocate Transient Moment Data

            Run(x).Response(j,k).Velocity  = Raw(j).Velocity( Run(x).Response(j,k).Idx );  % Allocate Transient Velocity Data
            Run(x).Response(j,k).Load      = Raw(j).Load( Run(x).Response(j,k).Idx );      % Allocate Transient Normal Force Data
            Run(x).Response(j,k).Pressure  = Raw(j).Pressure( Run(x).Response(j,k).Idx );  % Allocate Transient Pressure Data

%% Fitting Function
            [Run(x).Fit(j,k).FyResponseFit, ~, RelaxationLength] = ...
                    StepSteerFyResponseFit( Run(x).Response(j,k).Time - Run(x).Response(j,k).Time(1), ...
                    Run(x).Response(j,k).Force, Run(x).Response(j,k).Velocity );     

            %[MzResponseFit, ~] = StepSteerMzResponseFit(~);
            
            Run(x).Fit(j,k).RelaxationLength   = RelaxationLength;
            Run(x).Fit(j,k).Pressure           = Mesh(PressureCase).Pressure;
            Run(x).Fit(j,k).Load               = Mesh(j).Load;
            Run(x).Fit(j,k).Slip               = Mesh(1).Slip.Angle(SlipCase);
            
            SlipCase = SlipCase+1;
            if SlipCase > 3
                SlipCase = 1;
            end
        end
        
        PressureCase = PressureCase+1;
        if PressureCase > 3
            PressureCase = 1;
        end    
    end
end

%% Plotting Function
RelaxationLengthPlotting( Tire, Raw, Run, Idx, Mesh, Figure );