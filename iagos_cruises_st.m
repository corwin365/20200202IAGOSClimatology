clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%where is the data?
Settings.DataDir = [LocalDataDir,'/IAGOS/Timeseries'];
Settings.OutDir  = [LocalDataDir,'/corwin/IAGOS_ST/'];

%what time step to interpolate to?
Settings.dt = 5; %second

%how is a 'cruise' defined?
%first variable is height change required
%second variable is how long this change is over, in the time units specified above
Settings.MaxDz  = 100;%m
Settings.Window = 15*60./Settings.dt;  %fifteen minutes

%detrending window size
Settings.Detrend = 30*60./Settings.dt;

%smoothing size
Settings.Smooth = 5*60./Settings.dt;

%dates to loop over. A separate file will be produced for each day.
Settings.TimeScale = datenum(1994,1,1):1:datenum(2020,12,31);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% variables to retain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%geolocation
Settings.Vars.Geo.In = {'lat','lon','Time'};
Settings.Vars.Geo.Out  = {'Lat','Lon','Time'};

%s-transform
Settings.Vars.ST.In  = {'A','F1'};
Settings.Vars.ST.Out = {'A','k'};

%metadata
Settings.Vars.Meta.In  = {'air_press_AC','gps_alt_AC','zon_wind_AC','mer_wind_AC','air_temp_PM'};
Settings.Vars.Meta.Out = {'Prs','Z','U','V','T'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Settings.Window = round(Settings.Window);
if ~isodd(Settings.Window); Settings.Window = Settings.Window+1; end

Settings.Detrend = round(Settings.Detrend);
if ~isodd(Settings.Detrend); Settings.Detrend = Settings.Detrend+1; end

Settings.Smooth = round(Settings.Smooth);
if ~isodd(Settings.Smooth); Settings.Smooth = Settings.Smooth+1; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iDay = 1:1:numel(Settings.TimeScale)
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %create results arrays for day
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  X = NaN(1,1);%will expand out as needed
  Results = struct();
  Vars = [Settings.Vars.Geo.Out,Settings.Vars.ST.Out,Settings.Vars.Meta.Out];%,Settings.Vars.Special];
  for iVar=1:1:numel(Vars);
    Results.(Vars{iVar}) = X;
  end
  clear iVar Vars X
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %load day's data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  [yy,mm,dd] = datevec(Settings.TimeScale(iDay));
  disp(datestr(Settings.TimeScale(iDay)))
  
  
  
  Files = wildcardsearch([Settings.DataDir,...
                          '/',sprintf('%02d',yy), ...
                          '/',sprintf('%02d',mm),'/'],...
                          ['*',sprintf('%04d',yy),sprintf('%02d',mm),sprintf('%02d',dd),'*.nc']);
  for iFile=1:1:numel(Files);
    
    
    try
      %load file, including unit conversions
      Data = prep_iagos(Files{iFile}, ...
                        'SamplingRate',1./24./60./60.*Settings.dt, ...
                        'CruiseDz',Settings.MaxDz, 'CruiseWindow',Settings.Window, ...
                        'ApplyFlags',true);
      
      %loop over cruises
      for iCruise = 1:1:size(Data.Cruises,1)
        
        %identify cruise
        Cruise = Data.Cruises(iCruise,:);
        Cruise = Cruise(~isnan(Cruise));
        
        %extract data
        Var = Data.air_temp_AC(Cruise);
        Lon = Data.lon(Cruise);
        Lat = Data.lat(Cruise);
        
        %smooth out small scales
        Var = smooth(Var,Settings.Smooth);
        
        %detrend large scales
        BG = smooth(Var,Settings.Detrend);
        Var = Var - BG;
        
        %s-transform
        ST = nph_ndst(Var);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %variable storage
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for VarType = {'Geo','ST','Meta'}
          
          for iVar=1:1:numel(Settings.Vars.(VarType{1}).Out)
            
            %pull results array out of struct
            R = Results.(Settings.Vars.(VarType{1}).Out{iVar});
            
            %get data to store
            switch VarType{1};
              case 'ST';  O = ST;   C = 1:numel(Cruise);
              otherwise;  O = Data; C = Cruise;
            end
            O = O.(Settings.Vars.(VarType{1}).In{iVar});
            
            %store data, and fill any end-gaps with NaNs
            oldsize = size(R,2);
            R(iCruise,1:numel(C)) = O(C);
            if oldsize > numel(C);
              R(iCruise,numel(C)+1:oldsize) = NaN;
            end
            
            %put results back into struct
            Results.(Settings.Vars.(VarType{1}).Out{iVar}) = R;
          end; clear R O iVar oldsize
          
        end
      end; clear iFile
      
      
      
    catch; end
  end; clear iFile
  
  %finally, store the data for the day
  OutFile = [Settings.OutDir,'/IAGOS_ST_',num2str(Settings.TimeScale(iDay)),'.mat'];
  save(OutFile,'Results','Settings')

end