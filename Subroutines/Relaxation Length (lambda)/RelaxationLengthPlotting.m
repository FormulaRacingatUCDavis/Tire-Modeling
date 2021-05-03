function RelaxationLengthPlotting( Tire, Raw, Response, Idx, Mesh, Figure )
% Plots results from Relaxation Length (RL) Fitting

%% Plot the cases of the fit
PlotNumber=1;
for j = 1 : numel(Raw)
    for k = 1 : numel(Idx.NewBreaks)
        [FyResponseFit, ~] = ...
            StepSteerFyResponseFit( Response(j,k).Time - Response(j,k).Time(1), ...
            Response(j,k).Force );
        %{
        [MzResponseFit, ~] = ...
            StepSteerMzResponseFit( ( Response(j,k).Time - Response(j,k).Time(1) ), ...
            Response(j,k).Moment );
        %}
            
        figure(1);
        
        ax1 = nexttile;
        title('RL Fy Plot #' + string(PlotNumber));
        
        yyaxis(ax1, 'left');
        scatter( Response(j,k).Time - Response(j,k).Time(1), ...
            Response(j,k).Force );
        hold on;
        plot( FyResponseFit{2,1} );
        xlabel( 'Time: $t$' );
        ylabel( 'Cornering Force: $F_{y}$' );
        
        ax1.XAxisLocation = 'origin';
        if abs( max( Response(j,k).Force ) ) > abs( min( Response(j,k).Force ) )
            ylim( [ -abs( max( Response(j,k).Force ) ) abs( max( Response(j,k).Force ) ) ] );
        elseif abs( max( Response(j,k).Force ) ) < abs( min( Response(j,k).Force ) )
            ylim( [ -abs( min( Response(j,k).Force ) ) abs( min( Response(j,k).Force ) ) ] );
        else
            ylim auto;
        end
        
        figure(2);
        
        ax2 = nexttile;
        title('RL Mz Plot #' + string(PlotNumber));
        
        yyaxis(ax2, 'left');
        scatter( Response(j,k).Time - Response(j,k).Time(1), ...
            Response(j,k).Moment );
        %hold on;
        %plot( MzResponseFit{2,1} );
        xlabel( 'Time: $t$' );
        ylabel( 'Aligning Moment: $M_{z}$' );
        
        ax2.XAxisLocation = 'origin';
        if abs( max( Response(j,k).Moment ) ) > abs( min( Response(j,k).Moment ) )
            ylim( [ -abs( max( Response(j,k).Moment ) ) abs( max( Response(j,k).Moment ) ) ] );
        elseif abs( max( Response(j,k).Moment ) ) < abs( min( Response(j,k).Moment ) )
            ylim( [ -abs( min( Response(j,k).Moment ) ) abs( min( Response(j,k).Moment ) ) ] );
        else
            ylim auto;
        end

        PlotNumber = PlotNumber+1;
    end
end

%% Surface plot Relaxation Length as a function of Pressure and Load
a=1;