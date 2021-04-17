function RelaxationLengthPlotting( Tire, Raw, Response, Idx, Mesh, Figure )
% Plots results from Relaxation Length (RL) Fitting

%% Plot the cases of the fit
PlotNumber=1;
for j = 1 : numel(Raw)
    for k = 1 : numel(Idx.NewBreaks)
        [FyResponseFit, ~] = ...
            StepSteerFyResponseFit( Response(j,k).Time, Response(j,k).Force );
        
        subplot( 6, 5, PlotNumber );
        xlabel( 'Time: $t$' );

        yyaxis left;
        scatter( Response(j,k).Time, Response(j,k).Force );
        hold on;
        plot( FyResponseFit{1} );
        ylabel( 'Cornering Force: $F_{y}$' );
        
        yyaxis right;
        scatter( Response(j,k).Time, Response(j,k).Moment );
        xlabel( 'Time: $t$' );
        ylabel( 'Aligning Moment: $M_{z}$' );
        
        PlotNumber = PlotNumber+1;
    end
end

%% Surface plot Relaxation Length as a function of Pressure and Load
PlotNumber=1;