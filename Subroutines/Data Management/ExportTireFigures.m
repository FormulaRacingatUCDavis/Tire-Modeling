function ExportTireFigures( Tire, Data, Bin, Figure, Directory )
%% Clear / Create Main Export Directory
ExportPath = [Directory.Media, '\', Tire.Name, ' (', Tire.Date, ')' ];

if exist( ExportPath, 'dir' )
    rmdir( ExportPath, 's' )
end

mkdir( ExportPath )

%% Export Developer Figures
mkdir( [ExportPath, '\Developer Figures'] )
StructureFigureSearch( Figure )

%% Pacejka Carpet Plots
%%% Slip & Load Carpet Plots
a = 1;

%%% Camber & Load Carpet Plots

%% Radial Deflection Carpet Plots

%% Local Functions
    function StructureFigureSearch( Struct )
        Field = fieldnames( Struct );
        
        for i = 1:numel( Field )
           switch class( Struct.(Field{i}) )
               case 'struct'
                   StructureFigureSearch( Struct.(Field{i}) )
               case 'matlab.ui.Figure'
                   ExportDeveloperFigure( Struct.(Field{i}) )
           end
        end
    end

    function ExportDeveloperFigure( Figure )
        for i = 1 : numel( Figure )
            Figure(i).WindowState = 'maximize';
            saveas( Figure(i), [ExportPath, '\Developer Figures\', Figure(i).Name], 'png' );
            Figure(i).WindowState = 'minimize';
        end
    end
end
    