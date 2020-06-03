clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% take AMDAR data from the US MDCRS programme, and reformat it as 
% IAGOS data to allow us to use common code for IAGOS analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Settings.DataDir   = [LocalDataDir,'/MDCRS/'];
Settings.OutDir    = [LocalDataDir,'/IAGOS/AMDAR_as_IAGOS/'];
Settings.TimeScale = datenum(2020,1,1):1:datenum(2020,1,16);


Settings.InVars = {'latitude','longitude','altitude','timeObs','temperature','windDir','windSpeed','ID'}; %'ID' is necessary for the programmatic logic

%the data is of dreadful quality. So we need to filter it heavily. 
%A lot of data is going to be lost, but we'll cope...
Settings.MinFlightPoints  = 50; %discard short segments of data
Settings.Maxdt            = 30;% minutes
Settings.Maxdx            = 200; %km
Settings.MaxSpeed         = 1000; %km/h

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iDay=1:1:numel(Settings.TimeScale)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% find all files for this day, then load and merge them
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  [yy,mm,dd] = datevec(Settings.TimeScale(iDay));
  DayPath = [Settings.DataDir,'/',num2str(yy),'/', ...
             sprintf('%04d',yy),sprintf('%02d',mm),sprintf('%02d',dd)];
  Files = wildcardsearch(DayPath,'*');
  clear DayPath
  clear Store
  for iFile=1:1:numel(Files)
    %load the file
    try
      FileData = rCDF(Files{iFile});
    catch
      continue
    end
    
    %produce a unique identifier for each flight
    %these may span multiple files, so we need to choose something
    %that can be reproduced between files
    TNU = char(FileData.en_tailNumber');
    ORI = char(FileData.origAirport');
    DES = char(FileData.destAirport');
    ORI = ORI(:,1:3); DES =  DES(:,1:3);
    ID  = NaN(size(TNU,1),1);
    for iID = 1:1:numel(ID);
      a = 100.*sum(double(TNU(iID,:)));
      b =  10.*sum(double(ORI(iID,:)));
      c =   1.*sum(double(DES(iID,:)));
      ID(iID) = a+b+c;
    end
    FileData.ID = ID;
    clear TNU ORI DES iID a b c ID

    %extract the data
    if ~exist('Store','var')
      for iVar=1:1:numel(Settings.InVars)
        Store.(Settings.InVars{iVar}) = FileData.(Settings.InVars{iVar});
      end
    else
      Store = cat_struct(Store,FileData,1);
    end
    clear iVar
    
  end; clear iFile FileData yy mm dd Files
  
  %drop any points where all the science variables are NaN
  Good = find((~isnan(Store.temperature)+~isnan(Store.windSpeed)) > 0);
  Store = reduce_struct(Store,Good);
  clear Good
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% split the data from points into flights
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  FlightData = struct();
  Flights = unique(Store.ID);
  Count = 0;

  for iFlight=1:1:numel(Flights)
    
    %find flight
    ThisFlight = find(Store.ID == Flights(iFlight));
    
    %simple checks - too long or too short?
    if numel(ThisFlight) < Settings.MinFlightPoints;continue; end
    
    %flight now needs ordering in time
    [~,idx] = sort(Store.timeObs(ThisFlight),'asc');
    ThisFlight = ThisFlight(idx);
    
    %any large time discontinuities? This indicates either data gaps or landings.
    if max(diff(Store.timeObs(ThisFlight)))./60 > Settings.Maxdx;continue; end
    
    stop
    
    %any large space discontinuities? This indicates either data gaps or landings.
    x1 = [Store.latitude(ThisFlight),Store.longitude(ThisFlight)];
    x2 = circshift(x1,1,1);
    dx = nph_haversine(x1,x2); dx = dx(2:end);%km
    if max(dx) > Settings.Maxdx; continue; end
    
    %any stupidly speedy? This indicates multiple flights in this "record" or bad data.
    dt = diff(Store.timeObs(ThisFlight))./60./60; %hours
    speed = dx./dt;
    if max(speed) > Settings.MaxSpeed; continue; end
    
    %ok. this is an acceptable flight. keep it.
    Count = Count+1;
    ThisFlight = reduce_struct(Store,ThisFlight);
    FlightData.(['Flight',sprintf('%04d',Count)]) = ThisFlight;
    
  end
  clear dt dx Flights idx iFlight speed x1 x2 ThisFlight Store
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% make these flights look like IAGOS data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  [yy,mm,dd] = datevec(Settings.TimeScale(iDay));
  FilePath = [Settings.OutDir,'/mdcrs_',sprintf('%04d',yy),'/'];
  for iFile=1:1:Count
    
    %pull out the data
    Data = FlightData.(['Flight',sprintf('%04d',iFile)]);
    
    %rename the variables to match IAGOS
    Out = struct();
    for iField=1:1:numel(Settings.InVars)
      switch Settings.InVars{iField}
        case 'longitude';   Out.lon = Data.longitude;
        case 'latitude';    Out.lat = Data.latitude;
        case 'altitude';    Out.baro_alt_AC = Data.altitude;
        case 'timeObs';     Out.UTC_time = (datenum(1970,1,1,0,0,Data.timeObs)-floor(datenum(1970,1,1,0,0,Data.timeObs))).*24.*60.*60;
        case 'temperature'; Out.air_temp_AC = Data.temperature;
        case 'windDir';
          Out.zon_wind_AC = Data.windSpeed .* cosd(90-Data.windDir);
          Out.mer_wind_AC = Data.windSpeed .* sind(90-Data.windDir);
        case 'windSpeed'; %already handled above
        case 'ID'; Out.flight_ID = Data.ID;
        otherwise; disp('Variable not defined'); stop;
      end
    end
    clear iField
    
    %outfile name
    %same date convention as IAGOS, with thei ID string replaced with mine
    OutFile = [FilePath,'/MDCRS_timeseries_',...
               sprintf('%04d',yy),sprintf('%02d',mm),sprintf('%02d',dd), ...
               num2str(Out.flight_ID(1)),'.nc'];
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %% make the minimally-accurate netCDF file to the job
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    %meta headers
    MetaData.title = 'MDCRS data reformatted as IAGOS';
    MetaData.long_title = '';
    % MetaData.comment = Settings.Comment;
    MetaData.CreatedBy = 'Corwin Wright, Univ. Bath, UK';
    cjw_nc_create(OutFile,MetaData,1);
    clear MetaData

    %dimensions
    Dimensions = struct();
    Dimensions = cjw_nc_prepop_dim(Dimensions,'Point',1:1:numel(Out.lat));
    Dimensions(end).FullName = 'Point number';
    cjw_nc_makedims(OutFile,Dimensions);
    clear Dimensions
    
    %fill with data
    Data = struct;
    for iField=1:1:numel(Settings.InVars)
      switch Settings.InVars{iField}
        case 'longitude';   Name = 'lon'; Units = 'degrees';
        case 'latitude';    Name = 'lat'; Units = 'degrees';
        case 'altitude';    Name = 'baro_alt_AC'; Units = 'metres';
        case 'timeObs';     Name = 'UTC_time'; Units = ['seconds since ',datestr(Settings.TimeScale(iDay),'yyyy-mm-dd'),' 00:00:00'];
        case 'temperature'; Name = 'air_temp_AC'; Units = 'K';
        case 'ID';          Name = 'flight_ID'; Units = '';          
          %next two are special cases
        case 'windDir';     Name = 'zon_wind_AC'; Units = 'm/s';
        case 'windSpeed';   Name = 'mer_wind_AC'; Units = 'm/s';
        otherwise; disp('Variable not defined'); stop;
      end
      
      Data = cjw_nc_prepop_data(Data,Name);
      Data.(Name).Dims  = [1];
      Data.(Name).Units = Units;
      Data.(Name).Data  = Out.(Name);
    end; clear iField Name Units
    
    %write
    cjw_nc_writedata(OutFile,Data);
    clear Data OutFile Data
    

    
  
  end; clear iFile Count dd yy mm FilePath
  
  disp(['Done ',datestr(Settings.TimeScale(iDay))]);  
end; clear iDay