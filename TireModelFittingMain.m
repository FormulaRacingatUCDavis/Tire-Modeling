clc; clear; clear('global'); close all;

%% FRUCD Tire Model Fitting Tool Main Script
% This script fits a Pacejka contact patch tire model and radial deflection 
% models to Calspan TTC Data Sets and creates an instance of the FRUCDTire 
% class for use in vehicle models.
% 
% Use MATLAB, SI Unit Calspan TTC Data Sets
% 
% Authors: 
% - Blake Christierson (Sep. 2018 - Jun. 2021) 
% - Carlos Lopez       (Jan. 2019 -          )

%% Initialization
% Sets up the model structure and adds relevant directories for saving and 
% storing results.

%%% Global Variables
global Directory Figure

%%% Figure Interpreter
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%%% Default Directories
Directory.Tool  = fileparts( matlab.desktop.editor.getActiveFilename );
Directory.Data  = [Directory.Tool(1:max(strfind( Directory.Tool,'\' ))), 'Tire-Data'];
Directory.Save  = [Directory.Tool, '\Models'];
Directory.Media = [Directory.Tool, '\Media'];

addpath( genpath( Directory.Tool ) );
addpath( genpath( Directory.Data ) );

%%% Figure Structure
Figure.Mode  = 'Debug';
Figure.State = 'minimized';

%%% Initializing Tire Model
Tire = ModelParameterSetup;

%% Data Import
% Imports and bins the FSAE TTC test data into Data and Bin structures which 
% are utilized throughout the rest of the fitting process to select data from 
% desired operating conditions. 

TestName = {'Transient'                 , ...
            'Cornering 1'               , ...
            'Cornering 2'               , ...
            'Warmup'                    , ...
            'Drive, Brake, & Combined 1', ...
            'Drive, Brake, & Combined 2'};

%%% Change Current Directory to Default Data Directory
cd( Directory.Data )

%%% Import Run Files
for i = 1 : 6
    [FileName, PathName] = uigetfile( { '*.mat;*.dat' }, ...
         ['Please select the ', TestName{i}, ' Run Data File'] );  
    
    if ~ischar( FileName )
        continue
    else    
        Data(i) = DataImport( fullfile( PathName, FileName ), TestName{i} ); %#ok<SAGROW>
    end
end

%%% Empty Check
if exist( 'Data', 'var' ) == 0
    error( 'No Run Data Files Selected' );
end

clear TestName FileName PathName

%%% Binning Run Data
for i = find( ~cellfun( @isempty, { Data(:).Source } ) )
   Bin(i) = DataBinning( Data(i) ); %#ok<SAGROW>
end

clear i

%% Contact Patch Load Modeling
%%% Steady State, Pure Slip Force & Aligning Moment Fitting
Tire = PureLongitudinalFitting( Tire, Data, Bin ); % Longitudinal Force ( Fxo )

Tire = PureLateralFitting( Tire, Data, Bin ); % Lateral Force ( Fyo )

% load('TestData_1_27_21_17_15.mat')

%Tire = PureAligningFitting( Tire, Data, Bin ); % Aligning Moment ( Mzo )

%%% Steady State, Combined Slip Force Modeling
% This is currently undeveloped due to data limitations. Instead, combined tire forces
% can be evaluated using the Modified-Nicolas-Comstock (MNC) Model on the pure slip 
% models. It is implemented within the FRUCDTire Class Definition.
Model = FRUCDTire([], [], [], Tire.Pacejka, [], [], []);

%%% Steady State, Combined Slip Moment Fitting
% Tire = CombinedAligningFitting( Tire, Data, Bin ); % Aligning Moment (Mz)

% Tire = OverturningFitting( Tire, Data, Bin ); % Overturning Moment (Mx)

% Tire = ResistanceModeling( Tire, Data, Bin ); % Rolling Resistance (My)

%%% Transient Response
% Tire = RelaxationLengthFitting( Tire, Data, Bin );

return

%% Radial Deflection Modeling
% Cornering Stiffness Modeling

%% Thermal Modeling


%% Exporting Model
