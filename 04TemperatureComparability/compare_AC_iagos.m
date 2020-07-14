clearvars


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare AC temp (which we use) to IAGOS package temp (which is better)
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/05/19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%file handling
Settings.DataDir = [LocalDataDir,'/IAGOS/Timeseries'];
Settings.TimeScale = datenum(2000,1,1):1:datenum(2000,1,31);

%variables we are comparing
Settings.Vars = {'air_temp_AC','air_temp_PM','air_stag_temp_AC','air_stag_temp_PM'};

%cruise identification
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%what time step to interpolate to?
Settings.dt = 5; %second

%how is a 'cruise' defined?
%first variable is height change required
%second variable is how long this change is over, in the time units specified above
Settings.MaxDz  = 100;%m
Settings.Window = 15*60./Settings.dt;  %fifteen minutes

%minimum length (km) of cruise
Settings.MinCruiseLength = 1000;

%spectral analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%distance spacing
Settings.SA.dx = 1; %km

%maximum gap size
Settings.SA.MaxSpaceGap = 40; %km

%low-pass filter size
Settings.SA.Detrend = 1000./Settings.SA.dx;

%high-pass filter size
Settings.SA.Smooth = 20./Settings.SA.dx;

%gridding for the comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%

Settings.Grid.TPrime = -5:0.05:5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create results array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Results = zeros(numel(Settings.Grid.TPrime),numel(Settings.Grid.TPrime),numel(Settings.Vars)-1);

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
  %load day's data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  [yy,mm,dd] = datevec(Settings.TimeScale(iDay));
  disp(datestr(Settings.TimeScale(iDay)))
  
  Files = wildcardsearch([Settings.DataDir,'/iagos_',sprintf('%04d',yy)],...
                         ['*',sprintf('%04d',yy),sprintf('%02d',mm),sprintf('%02d',dd),'*.nc']);
  clear yy mm dd
  for iFile=1:1:numel(Files);
    
    
    %load file, including unit conversions
    
    %in this step we interpolate to time to identify the cruises. space
    %intepolation will be done instead in the loop below to get accurate
    %spatial wavelengths
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
      dxS = dx; for iX=2:1:numel(dxS); dxS(iX) = dxS(iX-1)+dxS(iX); end; clear iX
      dx2 = 0:Settings.SA.dx:max(dxS);
      
      %check all distance elements are unique (drop any duplicates below)
      if find(dx == 0); [~,Unique] = unique(dxS,'stable');
      else              Unique = 1:1:numel(dxS);
      end
      
      %check we still have enough data to be useful
      if numel(Unique) < 10; continue; end
      
      %check the cruise is long enoigh
      if max(dxS) < Settings.MinCruiseLength; continue; end
      
      %interpolate all variables to fixed grid
      Vars = fieldnames(Data);
      Regular = struct();
      for iVar=1:1:numel(Vars)
        if strcmp(Vars{iVar},'MetaData');     continue; end
        if strcmp(Vars{iVar},'OriginalTime'); continue; end
        if strcmp(Vars{iVar},'Cruises');      continue; end
        
        Var = Data.(Vars{iVar});
        if sum(~isnan(Var(Cruise))) <2 ; continue; end
        Regular.(Vars{iVar}) = interp1gap(dxS(Unique),Var(Cruise(Unique)),dx2,Settings.SA.MaxSpaceGap);
        
      end
      Regular.dxS = dx2;
      clear iVar Vars dx dxS dx2 x y Lat Lon 
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %detrend each temperature variable, then store for comparison
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      TStore = NaN(numel(Settings.Vars),1); % the 1 will expand
      
      for iVar = 1:1:numel(Settings.Vars)
        
        if ~isfield(Regular,Settings.Vars{iVar}); continue; end
        
        
        T = Regular.(Settings.Vars{iVar});
        
        %smooth out small scales
        T = smooth(T,Settings.SA.Smooth);
        
        
        %detrend large scales
        T = T-smooth(T,Settings.SA.Detrend);
        
        %Store
        TStore(iVar,1:numel(T)) = T;
        
        clear T
        
      end; clear iVar
      TStore(TStore == 0) = NaN;
      
      %ok. grid each one up against AC
      for iLine=2:1:4;
        
        x = TStore(    1,:);
        y = TStore(iLine,:);
        [xi,yi] = meshgrid(Settings.Grid.TPrime,Settings.Grid.TPrime);
        Good = find(~isnan(x+y));
        Results(:,:,iLine-1) = Results(:,:,iLine-1) + bin2mat(x(Good),y(Good),ones(size(Good())),xi,yi,'@nansum');
        
      end
      
      

    end; clear iCruise
  end; clear iFile
  
  
  
end; clear iDay
save('comparison.mat','Settings','Results')
