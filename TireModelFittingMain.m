clc; clear; close all;

%% TireModelFittingMain - FRUCD Tire Modeling Main Script
% This script fits a Pacejka contact patch tire model, radial deflection 
% models, and estimates thermal parameters from Calspan TTC data sets and 
% creates an instance of the FRUCDTire class for use in vehicle models.
% 
% Inputs:
%   - SI Unit Calspan TTC Data Sets [via uigetfile()]
% 
% Outputs:
%   Directory - Data, Model, and Figure Path Locations
%   Data      - Parsed FSAE TTC Data
%   Bin       - Logical Binnings for Separating Operating Conditions
%   Tire      - Tire Model
%   Figure    - Stores Model Figures
%
% Author(s): 
% Blake Christierson (bechristierson@ucdavis.edu) [Sep 2018 - Jun 2021] 
% Carlos Lopez       (calopez@ucdavis.edu       ) [Jan 2019 -         ]
% 
% Last Updated: 02-May-2021

%% Initialization
% Sets up the model structure and adds relevant directories for saving and 
% storing results.

%%% Figure Interpreter
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%%% Default Directories
Directory.Tool       = fileparts( matlab.desktop.editor.getActiveFilename );
Directory.Data       = [Directory.Tool(1:max(strfind( Directory.Tool,'\' ))), 'Tire-Data'];
Directory.Resources  = [Directory.Tool(1:max(strfind( Directory.Tool,'\' ))), 'MATLAB-Resources'];

Directory.Model = [Directory.Tool, '\Vehicle Model Resources\Models'];
Directory.Media = [Directory.Tool, '\Media'];

addpath( genpath( Directory.Tool      ) );
addpath( genpath( Directory.Data      ) );
addpath( genpath( Directory.Resources ) );

%%% Figure Structure
Figure.Mode  = 'Debug';
Figure.State = 'minimized';

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

%% Initialize Tire Model
ModelName = inputdlg( sprintf('%s\n%s','Enter Tire Model Name in Following Format: ', ...
    '{Manufacturer} {Compound} {Diameter}x{Width}-{Rim Diameter}x{Rim Width}'), ...
    '', 1, {'Tire'} );

Tire = TireParameters( ModelName{1}, [Data.Source], []);

clear ModelName

%% Radial Deflection Modeling
Tire = RadialDeflectionFitting( Tire, Data ); % Radial Deflection (Re, Rl)

%% Contact Patch Load Modeling
%%% Steady State, Pure Slip Force & Aligning Moment Fitting
Tire = PureLongitudinalFitting( Tire, Data, Bin, Figure ); % Longitudinal Force ( Fxo )

Tire = PureLateralFitting( Tire, Data, Bin, Figure ); % Lateral Force ( Fyo )

Tire = PureAligningFitting( Tire, Data, Bin, Figure ); % Aligning Moment ( Mzo )

%%% Steady State, Combined Slip Force Modeling
% This is currently undeveloped due to data limitations. Instead, combined tire forces
% can be evaluated using the Modified-Nicolas-Comstock (MNC) Model on the pure slip 
% models. It is implemented within the FRUCDTire Class Definition.

%%% Steady State, Combined Slip Moment Fitting
% Tire = CombinedAligningFitting( Tire, Data, Bin, Figure ); % Aligning Moment (Mz)

Tire = OverturningFitting( Tire, Data, Bin, Figure ); % Overturning Moment (Mx)

% Tire = ResistanceModeling( Tire, Data, Bin, Figure ); % Rolling Resistance (My)

%%% Transient Response
% Tire = RelaxationLengthFitting( Tire, Data, Bin, Figure );

%% Thermal Modeling
% Heat Generation Modeling

%% Exporting Model
SaveModel = questdlg( 'Save Model?', '', 'Yes', 'No', 'No' );

if strcmpi( SaveModel, 'Yes' )
    save( [Directory.Model, '\', ...
        strrep(strrep(Tire.Name, '.', ''),' ', '_')], 'Tire' )
    
    ExportTireFigures( Tire, Data, Bin, Figure, Directory );
end

clear SaveModel