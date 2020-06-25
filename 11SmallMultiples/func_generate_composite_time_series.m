function func_generate_time_series(FILENAME, LON, LAT, RANGE, PRSRANGE,TIMESCALE, VARIABLE, METRIC, TIMEWINDOW)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate time series of GW properties measured by IAGOS
%
%Data are returned as an ANNUAL CYCLE
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/Jun/07
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
  R = NaN(366,...
          numel(Settings.Metrics));
  
  %extract the variable
  Var = Store.(Settings.Vars{iVar});
  
  %grid the data
%   textprogressbar(['Gridding ',Settings.Vars{iVar},' '])
  dd = date2doy(Store.Time);
  for iDay =1:1:366

    TimeWindow = iDay+[-1,1].*Settings.TimeWindow./2;
    
    InTimeRange= find(dd >= TimeWindow(1) & dd < TimeWindow(2));
    if min(TimeWindow) <   0; InTimeRange = [InTimeRange;find(dd > 366+min(TimeWindow))]; end
    if max(TimeWindow) > 366; InTimeRange = [InTimeRange;find(dd < max(TimeWindow)-366)]; end
    
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
      
%     textprogressbar(iTime./numel(Settings.TimeScale).*100)
  end
  
  %store
  Results.(Settings.Vars{iVar}) = R;
                                      
%   textprogressbar(100); textprogressbar('!');
  
  clear R Var iTime InTimeRange iMetric
  
end; clear iVar
save(Settings.OutFile,'Settings','Results')