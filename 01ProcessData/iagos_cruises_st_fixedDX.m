clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%general setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%where is the data?
Settings.DataDir = [LocalDataDir,'/IAGOS/Timeseries'];
Settings.OutDir  = [LocalDataDir,'/corwin/IAGOS_st/'];

%dates to loop over. A separate file will be produced for each day.
Settings.TimeScale = datenum(2018,1,1):1:datenum(2018,12,31);

%cruise identification
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%what time step to interpolate to?
Settings.dt = 5; %second

%how is a 'cruise' defined?
%first variable is height change required
%second variable is how long this change is over, in the time units specified above
Settings.MaxDz  = 100;%m
Settings.Window = 15*60./Settings.dt;  %fifteen minutes

%spectral analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%distance spacing
Settings.SA.dx = 1; %km

%low-pass filter size
Settings.SA.Detrend = 1000./Settings.SA.dx;

%high-pass filter size
Settings.SA.Smooth = 5./Settings.SA.dx;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% variables to retain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%geolocation
Settings.Vars.Geo.In  = {'lat','lon','Time'};
Settings.Vars.Geo.Out = {'Lat','Lon','Time'};

%s-transform - temperature
Settings.Vars.STT.In  = {'IN','A','F1'};
Settings.Vars.STT.Out = {'Tprime','STT_A','STT_k'};

%s-transform - U
Settings.Vars.STU.In  = {'IN','A','F1'};
Settings.Vars.STU.Out = {'Uprime','STU_A','STU_k'};

%s-transform - V
Settings.Vars.STV.In  = {'IN','A','F1'};
Settings.Vars.STV.Out = {'Vprime','STV_A','STV_k'};

%metadata
Settings.Vars.Meta.In  = {'air_press_AC','baro_alt_AC','zon_wind_AC','mer_wind_AC','air_temp_PM'};
Settings.Vars.Meta.Out = {'Prs','Z','U','V','T'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Settings.Window = round(Settings.Window);
if ~isodd(Settings.Window); Settings.Window = Settings.Window+1; end

Settings.SA.Detrend = round(Settings.SA.Detrend);
if ~isodd(Settings.SA.Detrend); Settings.SA.Detrend = Settings.SA.Detrend+1; end

Settings.SA.Smooth = round(Settings.SA.Smooth);
if ~isodd(Settings.SA.Smooth); Settings.SA.Smooth = Settings.SA.Smooth+1; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iDay = 1:1:numel(Settings.TimeScale)
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %create results arrays for day
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  X = NaN(1,1);%will expand out as needed
  Results = struct();
  Vars = [Settings.Vars.Geo.Out, ...
          Settings.Vars.STT.Out, ...
          Settings.Vars.STU.Out, ...
          Settings.Vars.STV.Out, ...          
          Settings.Vars.Meta.Out];
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
    
    
%    try
      %load file, including unit conversions
      Data = prep_iagos(Files{iFile}, ...
                        'SamplingRate',1./24./60./60.*Settings.dt, ...
                        'CruiseDz',Settings.MaxDz, 'CruiseWindow',Settings.Window, ...
                        'ApplyFlags',true);
                      
      %loop over cruises
      for iCruise = 1:1:size(Data.Cruises,1)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %find cruise info and regularly grid
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        %identify cruise
        Cruise = Data.Cruises(iCruise,:);
        Cruise = Cruise(~isnan(Cruise));
                
        %extract geolocation
        Lat = Data.lat(Cruise);    
        Lon = Data.lon(Cruise);     
        
        %get distances for regular interpolation
        x = [Lat;Lon]'; y = circshift(x,1,1);
        dx = nph_haversine(x,y); dx(1) = dx(2);
        dxS = dx; for iX=2:1:numel(dxS); dxS(iX) = dxS(iX-1)+dxS(iX); end  
        dx2 = 0:Settings.SA.dx:max(dxS);

        %check we still have enough data to be useful
        if numel(dx2) < 10; continue; end
        
        %interpolate all variables to fixed grid
        Vars = fieldnames(Data);
        Regular = struct();
        for iVar=1:1:numel(Vars)
          if strcmp(Vars{iVar},'MetaData');     continue; end
          if strcmp(Vars{iVar},'OriginalTime'); continue; end
          if strcmp(Vars{iVar},'Cruises');      continue; end
          
          Var = Data.(Vars{iVar});
          Regular.(Vars{iVar}) = interp1(dxS,Var(Cruise),dx2);
        end
        clear iVar Vars dx dxS dx2 x y
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %detrend and s-transform T, U and V
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        T = Regular.air_temp_AC;
        U = Regular.zon_wind_AC;
        V = Regular.mer_wind_AC;
        
        %smooth out small scales
        T = smooth(T,Settings.SA.Smooth);
        U = smooth(U,Settings.SA.Smooth);
        V = smooth(V,Settings.SA.Smooth);        
        
        %detrend large scales
        T = T-smooth(T,Settings.SA.Detrend);
        U = U-smooth(U,Settings.SA.Detrend);
        V = V-smooth(V,Settings.SA.Detrend);

        %s-transform
        STT = nph_ndst(T);
        STU = nph_ndst(U);
        STV = nph_ndst(V);
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %store the outputs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for VarType = {'Geo','STT','STU','STV','Meta'}
          
          for iVar=1:1:numel(Settings.Vars.(VarType{1}).Out)
           
            %pull results array out of struct
            R = Results.(Settings.Vars.(VarType{1}).Out{iVar});
            
            %get data to store
            switch VarType{1};
              case 'STT'; O = STT;
              case 'STU'; O = STU;
              case 'STV'; O = STV;
              otherwise;  O = Regular;
            end
            O = O.(Settings.Vars.(VarType{1}).In{iVar});
            
            %store data, and fill any end-gaps with NaNs
            oldsize = size(R,2);
            R(iCruise,1:numel(O)) = O;
            if oldsize > numel(O);
              R(iCruise,numel(O)+1:oldsize) = NaN;
            end
            
            %put results back into struct
            Results.(Settings.Vars.(VarType{1}).Out{iVar}) = R;
          end; clear R O iVar oldsize
          
        end
      end; clear iFile
      
      
      
%    catch; end
  end; clear iFile
  
  %finally, store the data for the day
  OutFile = [Settings.OutDir,'/IAGOS_ST_',num2str(Settings.TimeScale(iDay)),'_v3.mat'];
  save(OutFile,'Results','Settings')

end
