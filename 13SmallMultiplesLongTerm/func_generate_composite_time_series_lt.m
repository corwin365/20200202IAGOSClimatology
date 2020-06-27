function func_generate_time_series_lt(FILENAME, LON, LAT, RANGE, PRSRANGE,TIMESCALE, VARIABLE, METRIC, TIMEWINDOW)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate time series of GW properties measured by IAGOS
%
%Data are returned as an SINGLE TIME SERIES
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/Jun/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Settings.DataDir    = [LocalDataDir ,'/corwin/IAGOS_annual/'];
Settings.TimeScale  = TIMESCALE;
Settings.TimeWindow = TIMEWINDOW; %days sliding window

Settings.Vars = {VARIABLE};
Settings.Metrics = {METRIC};
Settings.MinLambda = 25; %km

Settings.Range = RANGE;
Settings.PrsRange = PRSRANGE;
Settings.OutFile = ['data/',FILENAME,'.mat'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load all data for the region into one pile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Store.OldYear = 0;

for Time = floor(min(Settings.TimeScale)-(Settings.TimeWindow./2)) ...
           :1: ...
           ceil(max(Settings.TimeScale)+(Settings.TimeWindow./2)); 

  [yy,~,~] = datevec(Time);
  if Store.OldYear ~= yy

    %load data
    FileName = [Settings.DataDir,'/merged_',num2str(yy),'_sgolay900.mat'];
    if ~exist(FileName,'file'); continue; end
    File = load(FileName); File = File.Results;
    
    %create an array for unique flight IDs
    File.NU = repmat(1:1:size(File.Lat,1),size(File.Lat,2),1)';

    
    %extract those in our space region (just keep the whole time range, useful for e.g. smoothing)
    x = [File.Lat(:),File.Lon(:)];
    y = repmat([LAT,LON],size(x,1),1);
    dx = nph_haversine(x,y);
    InRange = find(dx < Settings.Range);
    
    InPrsRange = inrange(File.Prs,Settings.PrsRange);
    InRange = intersect(InRange,InPrsRange);
    InLambdaRange = find(1./File.STT_k > Settings.MinLambda);
    InRange = intersect(InRange,InLambdaRange);
    
    InLambdaRange = find(File.U > 20);
    InRange = intersect(InRange,InLambdaRange);    
    
    %glue to our arrays
    if ~isfield(Store,'Lat'); Store = reduce_struct(File,InRange); 
    else;                     Store = cat_struct(Store,reduce_struct(File,InRange),1,{'OldYear'});
    end
    Store.OldYear = yy;
    clear InLatRange InLonRange InRange File InPrsRange InLambdaRange
    disp(['Loaded ',num2str(yy)]);
  end
              
end; clear TimeRange yy FileName Time
disp('Data loaded')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for iVar=1:1:numel(Settings.Vars)
  
  %generate a results array
  R = NaN(numel(Settings.TimeScale),...
          numel(Settings.Metrics));
  
  %extract the variable
  Var = Store.(Settings.Vars{iVar});
  if numel(Var) == 0; Results.(Settings.Vars{iVar}) = R; continue; end %no data in the region. usually should only fire in tests...
  
  
  %grid the data
  textprogressbar(['Gridding ',Settings.Vars{iVar},' '])
  for iDay =1:1:numel(Settings.TimeScale)

    TimeWindow = Settings.TimeScale(iDay)+[-1,1].*Settings.TimeWindow./2;
    
    InTimeRange= find(Store.Time >= TimeWindow(1) & Store.Time < TimeWindow(2));
    if numel(InTimeRange) == 0; continue; end 
    
    NUnique = numel(~isnan(unique(Store.NU(InTimeRange))));
    
    %do calculation
    for iMetric = 1:1:numel(Settings.Metrics)
      
      %gini coeff only used for wave amplitudes
      if strcmp(Settings.Metrics{iMetric},'ginicoeff')
        if ~strcmp(Settings.Vars{iVar},'STT_A'); 
          continue; 
        end
      end
      
      %discard bad data
      Series = Var(InTimeRange);
      Good = find(~isnan(Series));
      if nansum(Series(Good)) == 0; continue; end
      
      if     strcmp(Settings.Metrics{iMetric},'prctile'); R(iDay,iMetric) = feval(Settings.Metrics{iMetric},Series(Good),Settings.PC(iMetric));
      elseif strcmp(Settings.Metrics{iMetric},'NU');      R(iDay,iMetric) = NUnique;
      elseif strcmp(Settings.Metrics{iMetric},'N');       R(iDay,iMetric) = numel(Good);
      else                                                R(iDay,iMetric) = feval(Settings.Metrics{iMetric},Series(Good));
      end
    end; clear iMetric Good Series
      
    textprogressbar(iDay./numel(Settings.TimeScale).*100)
  end
  
  %store
  Results.(Settings.Vars{iVar}) = R;
                                      
  textprogressbar(100); textprogressbar('!');
  
  clear R Var iTime InTimeRange iMetric
  
end; clear iVar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% also generate some surface time series over the same area and timescale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iVar=1:1:3;
  
  switch iVar
    case 1; InFile = [LocalDataDir,'/Miscellany/era5_10m_wind.mat'];         VarName = 'u10';
    case 2; InFile = [LocalDataDir,'/Miscellany/era5_surface_pressure.mat']; VarName = 'sp';
    case 3; InFile = [LocalDataDir,'/Miscellany/era5_total_precip.mat'];     VarName = 'tp';
  end
  
  %load data
  Data = load(InFile);
  Data = Data.(VarName);
  
  %find all points in range
  [xi,yi] = meshgrid(Data.lon,Data.lat);
  sz = size(xi);
  xi = xi(:); yi = yi(:);
  dx = nph_haversine([xi,yi],repmat([LON;LAT],1,numel(xi))');
  InRange = find(dx < Settings.Range);
  
  %make a single time series for the region
  TimeSeries = NaN(numel(Data.time),1);
  for iTime=1:1:numel(TimeSeries);
    Map = Data.(VarName);
    Map = Map(:,:,iTime);
    Map = Map(InRange);
    TimeSeries(iTime) = nanmean(Map(:));
  end
  
  %interpolate onto our amplitude time series
  %since both datasets are ~ monthly, this should be fine
  Results.(VarName) = interp1(Data.time,TimeSeries,TIMESCALE);
  

  
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% done
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


save(Settings.OutFile,'Settings','Results')