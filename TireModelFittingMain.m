clc; clear; clear('global'); close all;

% Use MATLAB, SI Unit Calspan TTC Data Sets

%% Initialization
% Global Variables
global Directory Figure

% Figure Interpreter
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

% Default Directories
Directory.Tool = mfilename('fullpath');
Directory.Tool = Directory.Tool(1 : strfind(Directory.Tool, mfilename()) - 2);
Directory.Data = [Directory.Tool, '\Tire Data']; % Change to Personal Directory Where Tire Data is Stored
Directory.Save = [Directory.Tool, '\Models']; % Change to Personal Directory Where Model & Documentation Should be Saved

addpath( genpath( Directory.Tool ) );

% Figure Structure
Figure.State = 'minimized';

% Initializing Tire Model
Tire = ModelParameterSetup;

%% Data Import
TestName = { 'Transient', 'Cornering 1', 'Cornering 2', ...
    'Warmup', 'Drive, Brake, & Combined 1', 'Drive, Brake, & Combined 2' };

% Change Current Directory to Default Data Directory
cd( Directory.Data )

% Import Run Files
for i = 1 : 6
    [FileName, PathName] = uigetfile( { '*.mat;*.dat' }, ...
         ['Please select the ', TestName{i}, ' Run Data File'] );  
    
    if ~ischar( FileName )
        continue
    else    
        Data(i) = DataImport( fullfile( PathName, FileName ), TestName{i} ); %#ok<SAGROW>
    end
end

% Empty Check
if exist( 'Data', 'var' ) == 0
    error( 'No Run Data Files Selected' );
end

clear TestName FileName PathName

%% Binning Run Data
for i = find( ~cellfun( @isempty, { Data(:).Source } ) )
   Bin(i) = DataBinning( Data(i) ); %#ok<SAGROW>
end

clear i

%% Steady State Pure Slip Fitting
 
% Lateral Force Fitting ( Fyo )
Tire = PureLateralFitting( Tire, Data, Bin );

% Longitudinal Force Fitting ( Fxo )
% Tire = PureLongitudinalFitting( Tire, Data, Bin );

% Aligning Moment Fitting ( Mzo )
% Tire = PureAligningFitting( Tire, Data, Bin );

return

%% Steady State Combined Slip Fitting
% This is currently undeveloped due to data limitations. Instead, combined
% tire forces can be evaluated using the Modified-Nicolas-Com

%% Past This is Undeveloped

%% Camber Evaluation Plots
for p = 1:length(Tire.run{2}.binval.P)
    
    Fig{p}.GAM = figure('Name','Optimal Camber from Tread Temperature Differential');
    
    for c = 2:length(Tire.run{2}.binval.IA)
        if p == 1
            i = 3;
        else 
            i = 2;
        end
        
        validIdx = Tire.run{i}.bin.P(:,p) & Tire.run{i}.bin.IA(:,c);
      
        if sum(validIdx) > 25
            
            subplot(3,1,c-1)
            hold on
            
            plot(Tire.run{i}.dat.ET(validIdx),...
                Tire.run{i}.dat.TSTC(validIdx) - ...
                Tire.run{i}.dat.TSTO(validIdx),...
                Tire.run{i}.dat.ET(validIdx),...
                Tire.run{i}.dat.TSTC(validIdx) - ...
                Tire.run{i}.dat.TSTI(validIdx));
            
            title({strcat('$P = ',num2str(Tire.run{i}.binval.P(:,p)), ...
                '$, $\gamma = ',num2str(Tire.run{i}.binval.IA(:,c)), ...
                '$')}, 'Interpreter','latex')
            
            xlabel('Elapsed Time (sec)','Interpreter','latex')
            ylabel('Temperature Differential ($^{\circ}$C)','Interpreter','latex')
            set(gca,'TickLabelInterpreter','latex')
            
            legend({'Outer Temperature Differential', ...
                'Inner Temperature Differential'},'Interpreter','latex')
            
            ylim([-15 15])
        end
    end
end

%% Radial Deflection Modeling
% Cornering Stiffness Modeling
for p = 2
    if p == 1
        i = 3;
    else
        i = 2;
    end
    
    [Tire.RL{1},Tire.RL{2}] = LR_P6( ...
        (16/2)*2.54 - ...
        Tire.run{i}.dat.RL(Tire.run{i}.bin.P(:,p)), ...
        abs(Tire.run{i}.dat.FZ(Tire.run{i}.bin.P(:,p))));
end

%% Transient Tire Dynamics
% Tire = TransientRelaxationLength( Tire, Data, Bin );