function Data = DataImport( FileID, TestName )
% Loads Calspan TTC Run Data Into Structure

[PathName , FileName , Ext ] = fileparts(FileID);

if strcmp( Ext, '.mat')   
    Raw = load(FileID);
    
    %%% Identification Data
    Data.Tire            = Raw.tireid;
    Data.TestName        = TestName;
    Data.Source.FileName = FileName;
    
    %%% String Parsing for Round and Run
    RoundIdx     = strfind( lower(PathName), 'round ');
    BackslashIdx = strfind( lower(PathName(RoundIdx:end)), '\'); 
    RunIdx       = strfind( lower(FileName), 'run' );

    Data.Source.Round = str2double( PathName( RoundIdx + (6 : BackslashIdx(1)-2) ) );
    Data.Source.Run   = str2double( FileName(RunIdx+3:end) );
    
    %%% Data Parsing
    % Note that FSAE TTC uses the SAE Z-Down coordinate system while FRUCD
    % used the SAE Z-Up coordinate system. The values labels as negative
    % are being converted between the coordinate systems.
    
    Data.Time = Raw.ET'; % ET: Elapsed Time [s]
    
    Data.Pressure    = Raw.P' * 0.145037737730217; % P: Pressure [psi]
    Data.Inclination = Raw.IA'; % IA: Inclination Angle [deg]
    Data.Velocity    = Raw.V';  % V: Road Velocity [m/s]
    Data.Omega       = Raw.N';  % N: Wheel Rotational Speed [rpm]
    
    Data.Slip.Angle = -Raw.SA'; % SA: Slip Angle [deg]
    Data.Slip.Ratio = Raw.SL';  % SL: Slip Ratio []
    
    Data.Force(1,:) = Raw.FX;  % FX: Longitudinal Force [N]
    Data.Force(2,:) = -Raw.FY; % FY: Lateral Force      [N] 
    Data.Force(3,:) = -Raw.FZ; % FZ: Normal Force       [N]
    
    Data.Mu(1,:) = Raw.NFX'; % NFX: Normalized Longitudinal Force ('Mu') []
    Data.Mu(2,:) = Raw.NFY'; % NFY: Normalized Lateral Force      ('Mu') []
    
    Data.Moment(1,:) = Raw.MX';  % MX: Overturning Moment [Nm]
    Data.Moment(3,:) = -Raw.MZ'; % MZ: Aligning Moment    [Nm]
    
    Data.Radius.Effective = Raw.RE'; % RE: Effective Radius [cm]
    Data.Radius.Loaded    = Raw.RL'; % RL: Loaded Radius    [cm]
     
    Data.Temp.Ambient = Raw.AMBTMP'; % AMBTMP: Ambient Temperature      [°C]
    Data.Temp.Surface = Raw.RST';    % RST   : Road Surface Temperature [°C]
    
    Data.Temp.Tire(1,:) = Raw.TSTI'; % TSTI: Inner Tire Surface Temperature  [°C]
    Data.Temp.Tire(2,:) = Raw.TSTC'; % TSTC: Center Tire Surface Temperature [°C]
    Data.Temp.Tire(3,:) = Raw.TSTO'; % TSTO: Outer Tire Surface Temperature  [°C]
    
    for ii = {'ET', 'P', 'IA', 'V', 'N', 'SA', 'SL', 'FX', 'NFX', 'MX', 'RE', 'AMBTMP'}
        if ~isfield( Data, 'Units') || isempty(Data.Units)
            Data.Units = Raw.channel.units(  strcmp( ii, Raw.channel.name) );
        else
            Data.Units(end+1) = Raw.channel.units(  strcmp( ii, Raw.channel.name) );
        end
    end
else
    error('Wrong File Type');
end

end
