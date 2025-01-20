function ExportTireFigures( Tire, Data, Bin, Figure, Directory )
%% Clear / Create Main Export Directory
ExportPath = [Directory.Media, '\', Tire.Name, ' (', Tire.Date, ')' ];

if exist( ExportPath, 'dir' )
    warning('off','all')
    rmdir( ExportPath, 's' )
    warning('on','all')
end

mkdir( ExportPath )

%% Export Developer Figures
mkdir( [ExportPath, '\Developer Figures'] )
StructureFigureSearch( Figure )

%% Pacejka Carpet Plots
%%% Slip Angle & Load Lateral Performance Carpet Plots
SplineFit = cell( numel( Bin(2).Values.Pressure )    , ...
                  numel( Bin(2).Values.Load )        , ...
                  numel( Bin(2).Values.Inclination ) );

PeakLine( numel( Bin(2).Values.Pressure )    , ...
                  numel( Bin(2).Values.Load )        , ...
                  numel( Bin(2).Values.Inclination ) ).Pacejka = []; 
for i = [2 3]
    for p = 1 : numel( Bin(2).Values.Pressure )       
        for c = 1 : numel( Bin(2).Values.Inclination )
            Idx.Valid = Bin(i).Pressure(p,:) & Bin(i).Inclination(c,:) & ...
                Bin(i).Slip.Ratio( find( Bin(i).Values.Slip.Ratio == 0 ), : );
            
            if all( sum( Idx.Valid & Bin(i).Load, 2 ) < 50 )
                    continue % Skip Sparse Bins
            end
            
            figure;
            
            % Scatter Raw Data
            RawScatter = scatter( rad2deg(Data(i).Slip.Angle( Idx.Valid )), Data(i).Force(2, Idx.Valid), 1, ...
                'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k',...
                'MarkerFaceAlpha',  .2, 'MarkerEdgeAlpha',  .2); 
            hold on;
            
            for z = 1 : numel( Bin(2).Values.Load ) 
                if sum( Idx.Valid & Bin(i).Load(z,:), 2 ) < 50 
                    continue % Skip Sparse Bins
                end
                
                % Plot Spline Fits
                SplineFit{p,z,c} = SmoothingSplineFit( ...
                    rad2deg(Data(i).Slip.Angle( Idx.Valid & Bin(i).Load(z,:) )), ...
                    Data(i).Force( 2, Idx.Valid & Bin(i).Load(z,:) ) );
                
                SplinePlot = plot( SplineFit{p,z,c}, 'b--' );
                
                % Plot Pacejka Fits
                SlipAngle = -15:0.01:15;
                
                Fidelity.Pure = 'Pacejka'; Fidelity.Combined = 'Pure';
                [~, FyPacejka, ~, ~, ~] = ContactPatchLoads( Tire, SlipAngle, ...
                    0, Bin(i).Values.Load(z), Bin(i).Values.Pressure(p), ...
                    Bin(i).Values.Inclination(c), 10, 1, ...
                    struct('Pure', 'Pacejka', 'Combined', 'MNC'));
                
                PacejkaPlot = plot( SlipAngle, FyPacejka, 'b', 'HandleVisibility', 'off' ); 
                
                % Label Data Set
                FySpline = feval( SplineFit{p,z,c}, SlipAngle )';
                % PlotInLineText( ['$F_{z}=', num2str(round(Bin(i).Values.Load(z))),'$ [$N$]'], ...
                %     [SlipAngle; FySpline], [14,1], 6 )
                
                % Find Line of Peaks
                [FyPacejkaMax, MaxIdx] = max( FyPacejka );
                [FyPacejkaMin, MinIdx] = min( FyPacejka );
                PeakLine(p,z,c).PacejkaSlip  = [SlipAngle( MaxIdx ); SlipAngle( MinIdx )];
                PeakLine(p,z,c).PacejkaForce = [FyPacejkaMax; FyPacejkaMin];
            end
            
            plot( [PeakLine(p,:,c).PacejkaSlip]', [PeakLine(p,:,c).PacejkaForce]', 'b' )
            
            title( { [Tire.Name, ' | Round ', num2str(Tire.Source(2).Round), ...
                ', Run ', num2str(Tire.Source(i).Run)], ...
                ['$P_{i}=', num2str( Bin(i).Values.Pressure(p) ), ...
                '$ [$psi$], $\gamma=', num2str( Bin(i).Values.Inclination(c) ), '$ [$deg$]'] } )
            
            xlabel( 'Slip Angle, $\alpha$ [$deg$]' ); xlim( [-17, 17] );
            ylabel( 'Lateral Force, $F_{y}$ [$deg$]' );
            
            legend( [RawScatter, SplinePlot, PacejkaPlot], ...
                {'Raw Data', 'Smoothing Spline', 'Pacejka MF6.2'} ) 
        end
    end
end

%%% Camber & Load Lateral Performance Carpet Plots

%%% Slip Ratio & Load Longitudinal Performance Carpet Plots

%% Radial Deflection Carpet Plots

%% Local Functions
    function StructureFigureSearch( Struct )
        Field = fieldnames( Struct );
        
        for ii = 1:numel( Field )
           switch class( Struct.(Field{ii}) )
               case 'struct'
                   StructureFigureSearch( Struct.(Field{ii}) )
               case 'matlab.ui.Figure'
                   ExportDeveloperFigure( Struct.(Field{ii}) )
           end
        end
    end

    function ExportDeveloperFigure( Figure )
        for ii = 1 : numel( Figure )
            Figure(ii).WindowState = 'maximize';
            saveas( Figure(ii), [ExportPath, '\Developer Figures\', Figure(ii).Name], 'png' );
            Figure(ii).WindowState = 'minimize';
        end
    end

    function Fit = SmoothingSplineFit( x, y )
        [xData, yData] = prepareCurveData( x, y );

        Form = fittype( 'Smoothingspline' );
        Opts = fitoptions( 'Method', 'SmoothingSpline' );
        Opts.Normalize = 'on';
        Opts.SmoothingParam = 0.95;

        Fit = fit( xData, yData, Form, Opts );
    end
end
    