function RelaxationLengthPlotting( Tire, Raw, Response, Idx, Mesh, Figure )
% Plots results from Relaxation Length (RL) Fitting

%% Plotting The Cases of The Fit
PlotNumber=1;
PressureCase=1;
SlipCase=1;

for j = 1 : numel(Raw)
    for k = 1 : numel(Idx.NewBreaks)
        [FyResponseFit, ~, RelaxationLength] = ...
            StepSteerFyResponseFit( Response(j,k).Time - Response(j,k).Time(1), ...
            Response(j,k).Force, Response(j,k).Velocity );     
        
        %[MzResponseFit, ~] = StepSteerMzResponseFit(~);
        
        figure(1);
        ax1 = nexttile;
                   
        title( { ['Pressure ($P$): $', num2str(Mesh(PressureCase).Pressure), '$ [$psi$]'], ...
             ['Load ($F_Z$): $', num2str(Mesh(j).Load), '$ [$lbf$]'], ...
             ['Slip Angle ($\alpha$): $', num2str(Mesh(1).Slip.Angle(SlipCase)), '$ [$rad$]'], ...
             ['Relaxation Length ($RL$): $', num2str(RelaxationLength), '$ [$m$]'] } );
         
        %{
        Using mode(Response(j,k).Slip) would provide the slip angle for 
        each cases of fit, but the reason why the the sign of the slip angle 
        is opposite to what was mentioned in the content round PDF is 
        because we use SAE z-up coordinate system, whereas FSAETTC uses
        ISO. From the graph, it would not show up because
        Mesh(1).Slip.Angle(SlipCase) contains the correct slip angle for
        each cases of fit. Hence, it is used instead.
        %}
        
        yyaxis(ax1, 'left');
        scatter( Response(j,k).Time - Response(j,k).Time(1), ...
            Response(j,k).Force );
        hold on;
        plot( FyResponseFit{1,1} );
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
        
        %{
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
        %}
        
        Fit.RelaxationLength(j,k)   = RelaxationLength;
        Fit.Pressure(j,k)           = Mesh(PressureCase).Pressure;
        Fit.Load(j,k)               = Mesh(j).Load;
        Fit.Slip(j,k)               = Mesh(1).Slip.Angle(SlipCase);       
        
        SlipCase = SlipCase+1;
        if SlipCase > 3
            SlipCase = 1;
        end
        
        PlotNumber = PlotNumber+1;
    end
    
    PressureCase = PressureCase+1;
    if PressureCase > 3
        PressureCase = 1;
    end
end

%% Surface Plot of Relaxation Length as a Function of Pressure and Load for Different Slip Angles
for i = 1:3
    figure(3);
    subplot(3,1,i);
     
    surf( ...
        reshape( Fit.Pressure(:,i), 3, 3 ), ...
        reshape( Fit.Load(:,i), 3, 3 ), ...
        reshape( Fit.RelaxationLength(:,i), 3, 3) ...
    );

    xlabel('Pressure: $P$'); 
    ylabel('Load: $F_{z}$'); 
    zlabel('Relaxation Length: $RL$');
    
    title( ['Slip Angle ($\alpha$): $', num2str( Fit.Slip(1,i) ), '$ [$rad$]'] );
end