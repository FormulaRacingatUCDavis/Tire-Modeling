clc; clear; clear('global'); close all;

%% FRUCD Tire Model Fitting Tool Main Script
% This script fits a Pacejka Tire Model to Calspan TTC Data Sets and
% creates an instance of the FRUCDTire class for use in vehicle models.
%
% Use MATLAB, SI Unit Calspan TTC Data Sets
%
% Authors: - Blake Christierson (Sept. 2018 - Jun. 2021)
%          - Carlos Lopez       (Jan.  2019 - ???      )

%% Initialization
% Global Variables
global Directory Figure

% Figure Interpreter
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

% Default Directories
Directory.Tool = fileparts( which( 'TireModelFittingMain.m' ) );
Directory.Data = [Directory.Tool(1:max(strfind( Directory.Tool,'\' ))), 'Tire-Data'];
Directory.Save = [Directory.Tool, '\Models'];

addpath( genpath( Directory.Tool ) );
addpath( genpath( Directory.Data ) );

% Figure Structure
Figure.State = 'minimized';

% Initializing Tire Model
Tire = ModelParameterSetup;

%% Data Import

%{
TestName = { 'Transient'                  , ...
             'Cornering 1'                , ...
             'Cornering 2'                , ...
             'Warmup'                     , ...
             'Drive, Brake, & Combined 1' , ...
             'Drive, Brake, & Combined 2' };

%%% Import Run Files
cd( Directory.Data )
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
if exist( 'Data', 'var' ) == 0; error( 'No Run Data Files Selected' ); end

%%% Binning Run Data
for i = find( ~cellfun( @isempty, { Data(:).Source } ) )
   Bin(i) = DataBinning( Data(i) ); %#ok<SAGROW>
end
%}

load('TestDataSet.mat');

%% Camber Evaluation Plots
figure;

subplot(2,2,1); hold on;
scatter( Data(5).Camber, Data(5).Mu(1,:), 3, 'b.', ...
    'MarkerFaceAlpha', 0.05, 'MarkerEdgeAlpha', 0.05 )

subplot(2,2,2); hold on;
scatter( Data(2).Camber, Data(2).Mu(2,:), 3, 'k.' )

subplot(2,2,3); hold on;
scatter( Data(2).Camber, Data(2).Temp.Tire(1,:)-Data(2).Temp.Tire(2,:) );

<<<<<<< Updated upstream
figure
ax1 = subplot(2,1,1); hold on;
scatter( Data(5).Temp.Tire(2,:), abs(Data(5).Mu(1,:)), 3, 'b.', ...
    'MarkerFaceAlpha', 0.05, 'MarkerEdgeAlpha', 0.01 )
bnd = boundary( Data(5).Temp.Tire(2,:)', abs(Data(5).Mu(1,:)'), 0.5 );
plot(  Data(5).Temp.Tire(2,bnd), abs(Data(5).Mu(1,bnd)), 'b' )
=======
return

Tire = PureLateralFitting( Tire, Data, Bin, Figure ); % Lateral Force ( Fyo )
>>>>>>> Stashed changes

ax2 = subplot(2,1,2); hold on;
scatter( Data(2).Temp.Tire(2,:), abs(Data(2).Mu(2,:)), 3, 'b.', ...
    'MarkerFaceAlpha', 0.05, 'MarkerEdgeAlpha', 0.01 )
bnd = boundary( Data(2).Temp.Tire(2,:)', abs(Data(2).Mu(2,:)'), 0.5 );
plot(  Data(2).Temp.Tire(2,bnd), abs(Data(2).Mu(2,bnd)), 'b' )

linkaxes( [ax1 ax2], 'x' )

return

%% Vertical (Radial) Stiffness
Tire = RadialDeformationFitting( Tire, Data );

%% Steady State, Pure Slip Contact Patch Force & Moments 
%%% Lateral Force Fitting ( Fyo )
Tire = PureLateralFitting( Tire, Data, Bin );

%%% Longitudinal Force Fitting ( Fxo )
Tire = PureLongitudinalFitting( Tire, Data, Bin );

%%% Aligning Moment Fitting ( Mzo )
Tire = PureAligningFitting( Tire, Data, Bin );

%% Steady State, Combined Slip Contact Patch Force & Moments
% This is currently undeveloped due to data limitations. Instead, combined
% tire forces can be evaluated using the Modified-Nicolas-Comstock (MNC)
% Model to couple pure slip model. It is implemented within the FRUCDTire 
% Class Definition.

%% Transient Tire Dynamics
% Tire = TransientRelaxationLength( Tire, Data, Bin );







%% Dt and Et Graphing
% Code does not affect fit of tire. Only Estimating Dt/Et parameters
dFz = @(Fz)(Fz - Tire.Fzo) ./ Tire.Fzo;
dPi = @(Pi)(Pi - Tire.Pio) ./ Tire.Pio;
figure
subplot(2,1,1)
gamma = 2;
Pi = 11;
Dt = @(Fz)(Tire.Ro .*(Fz./Tire.Fzo)) * (Tire.q.D.z(1) + Tire.q.D.z(2).*dFz(Fz) *...
    (1 - Tire.p.P.z(1) * dPi(Pi))) * (1 + Tire.q.D.z(3) * abs(gamma) + ...
    Tire.q.D.z(4) * gamma.^2);
 fplot(Dt,[0,2500])
 Ct = 0.8;
 Bt = 2.5;
 Slip = 0;
 subplot(2,1,2)
 Et = @(Fz)(Tire.q.E.z(1) + Tire.q.E.z(2)*dFz(Fz) + Tire.q.E.z(3) * dFz(Fz).^2)...
     *(1 + (Tire.q.E.z(4) + Tire.q.E.z(5)* gamma)*(2/pi)*atan(Bt * Ct * Slip))
 fplot(Et, [0,2500])