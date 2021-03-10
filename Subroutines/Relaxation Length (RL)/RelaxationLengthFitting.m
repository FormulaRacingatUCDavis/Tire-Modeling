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
% Last Updated: 10-Mar-2021

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
        Idx.Valid = Idx.Valid( Idx.Valid > find( Bin(1).Velocity(4,:), 1, 'last' ) );     % Finds last point of initial spring rate tests
        
        if numel( Idx.Valid ) < 50
            continue % Skip Sparse Bins
        end

        Raw(p,z).Slip       = Data(1).Slip.Angle(Idx.Valid);    % Allocate Slip Angle Data
        Raw(p,z).Force      = Data(1).Force(2,Idx.Valid);       % Allocate Lateral Force Data
        Raw(p,z).Velocity   = Data(1).Velocity(Idx.Valid);      % Allocate Velocity Data

        Raw(p,z).Pressure   = Data(1).Pressure(Idx.Valid);      % Allocate Pressure Data
        Raw(p,z).Load       = Data(1).Force(3,Idx.Valid);       % Allocate Normal Force Data

        Raw(p,z).dFz        = (Raw(p,z).Load     - Tire.Pacejka.Fzo) ./ Tire.Pacejka.Fzo;
        Raw(p,z).dPi        = (Raw(p,z).Pressure - Tire.Pacejka.Pio) ./ Tire.Pacejka.Pio;
    end
end

%% Filtering Data & Operating Conditions
Mesh(    ind2sub(size(Raw), find(cellfun(@isempty, {Raw.dFz}))) ) = [];
Raw(     ind2sub(size(Raw), find(cellfun(@isempty, {Raw.dFz}))) ) = [];

%% Velocity Thresholding, 
for j = 1: numel (Raw)
    Idx.IsRolling = find( Raw(j).Velocity > ( Bin(1).Tolerance.Velocity ./ 3 ) );   
end
%% Segmenting Velocity Based on Slope and Time Duration
for j = 1: numel (Rw)   
    %currently writing different for loops to divide steps
    
    %plot(diff(Idx.IsZero)) to identify where velocity is changing drastically
    %scatter(Idx.IsZero, Raw(j).Force(Idx.IsZero) ./ Raw(j).Load(Idx.IsZero)) )
    %to get a look at normalized force where this occur
    
    %Next steps:
    %Idx.IsPeaks = findpeaks();
    %Use NPeaks
    
    %{
    Idx.IsSteadyState = find( Threshold(p,z).Velocity > ...
        ( max( Threshold(p,z).Velocity ) - ( Bin(1).Tolerance.Velocity ./ 3 ) ) );
        
    SteadyState(p,z).Velocity = Threshold(p,z).Velocity(Idx.IsSteadyState);
    %}
    
    %{
    locating the beginning of change in velocities (represented by 
    first-order delays) through difference in the y-axis coordinate 
    (drastic change compared to the overall of the graph).
    %}
    
    %{
    IsSteadyState.Mean(p,z).Velocity = mean( Threshold(p,z).Velocity(Idx.IsSteadyState) );
    IsSteadyState.Min(p,z).Velocity = min( Threshold(p,z).Velocity(Idx.IsSteadyState) );
    IsSteadyState.DrasticChange(p,z).Velocity = IsSteadyState.Mean(p,z).Velocity - ...
        IsSteadyState.Min(p,z).Velocity;

    for i = 1: numel( Idx.IsSteadyState )
        slope = SteadyState(p,z).Velocity(i+1) - SteadyState(p,z).Velocity(i);
    
        if slope > IsSteadyState.DrasticChange.Velocity
            StartOfDrasticChangeInVelocity(i) = i; %error here
        end
    end
    %}
        
    %{
    time duration can be found by finding the difference between the
    previous beginning point and the next beginning point. Time
    duration that exceeds the average will be get rid of, leaving
    changes in velocities where transient-states occur
    %}
    
    %{
    for temp = 1: numel( StartOfDrasticChangeInVelocity )
    NumberOfDataPointsDuringSteadyState = ...
    StartOfDrasticChangeInVelocity(temp+1) - StartOfDrasticChangeInVelocity(temp);
    if NumberOfDataPointsDuringSteadyState < mean( NumberOfDataPointsDuringSteadyState )
    StartOfTransientStates(temp) = temp;
    end
    %}  
end
%% Segmenting Velocity for Out Transients
%{
%%%will be decided whether this step is necessary or not
for p = 1 : numel( Case.Pressure )       
    for z = 1 : numel( Case.Load )
        %{
        from the chosen points indicating beginning of change in velocities, outs of 
        the transient-states can be distinguished from the backs through 
        finding the point of maximum y-coordinate before having a negative 
        slope. (probably use diff() from one point to another to find the 
        slope and if negative, stop processing data and indicate 
        previously-calculated as "final" point of transient-state
        %}
                
        %{
        for more precise modeling, few iterations after the "final" point
        of transient-state will be included for plotting 
        %}
        
        %{
        keep in mind to include first velocity data points that are beyond
        defined threshold when modeling
        %}
        b=1;
    end
end
%}
%% Plotting Function
RelaxationLengthPlotting( Tire, Raw, Mesh, Figure );
