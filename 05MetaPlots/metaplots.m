clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%prepare IAGOS postprocessed metadata
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/05/14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%file handling
Settings.DataDir = [LocalDataDir,'/corwin/IAGOS_st'];
Settings.OutFile = 'metadata_all.mat';

%time period to include
Settings.TimeScale = datenum(1994,1,1):1:datenum(2020,12,31);


%map data
Settings.Grid.Lon = -180:5:180;
Settings.Grid.Lat = -90:5:90;
Settings.Grid.dTP = -200:5:200; %tropopause relative pressure (hPa)
Settings.Grid.Prs = 180:20:400;  %absolute pressure (hPa)

%time data
Settings.Grid.Years  = 1994:1:2020;
Settings.Grid.Months = 1:1:12; %yeah yeah

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prepare save grids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%maps
Results.Map.dTP  = zeros(numel(Settings.Grid.Lon),numel(Settings.Grid.Lat),numel(Settings.Grid.dTP));
Results.Map.Prs  = zeros(numel(Settings.Grid.Lon),numel(Settings.Grid.Lat),numel(Settings.Grid.Prs));

%time
Results.Time.dTP = zeros(numel(Settings.Grid.Years),numel(Settings.Grid.Months),numel(Settings.Grid.dTP));
Results.Time.Prs = zeros(numel(Settings.Grid.Years),numel(Settings.Grid.Months),numel(Settings.Grid.Prs));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get the data and store it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for iDay=1:1:numel(Settings.TimeScale)
  
  disp(datestr(Settings.TimeScale(iDay)))
  
  %get data
  File = [Settings.DataDir,'/IAGOS_ST_',num2str(Settings.TimeScale(iDay)),'_v3.mat'];
  if ~exist(File,'file'); clear File; continue; end
  Data = load(File); clear File
  
  %pull out data we need, and tidy it up
  Good = find(~isnan(Data.Results.STT_A) & Data.Results.STT_A~= 0);
  if numel(Good) == 0; clear Good Data; continue; end
  
  [Year,Month,~] = datevec(Data.Results.Time(Good));  
  Lon = Data.Results.Lon(Good);
  Lat = Data.Results.Lat(Good);
  TP  = Data.Results.TropPres(Good);
  Prs = Data.Results.Prs(Good);
  
  %produce TP-relative pressures
  dTP = Prs-TP; clear TP
  
  %produce maps
    %absolute pressure
  [xi,yi,zi] = meshgrid(Settings.Grid.Lon,Settings.Grid.Lat,Settings.Grid.Prs);
  Good = find(~isnan(Lon+Lat+Prs));
  Results.Map.Prs = Results.Map.Prs+ permute(bin2matN(3,Lon(Good),Lat(Good),Prs(Good),ones(size(Good)),xi,yi,zi,'@nansum'),[2,1,3]);
    %relative pressure
  [xi,yi,zi] = meshgrid(Settings.Grid.Lon,Settings.Grid.Lat,Settings.Grid.dTP);
  Good = find(~isnan(Lon+Lat+dTP));
  Results.Map.dTP = Results.Map.dTP+ permute(bin2matN(3,Lon(Good),Lat(Good),dTP(Good),ones(size(Good)),xi,yi,zi,'@nansum'),[2,1,3]);
  
  %produce time grids
    %absolute pressure
  [xi,yi,zi] = meshgrid(Settings.Grid.Years,Settings.Grid.Months,Settings.Grid.Prs);
  Good = find(~isnan(Year+Month+Prs));
  Results.Time.Prs = Results.Time.Prs+ permute(bin2matN(3,Year(Good),Month(Good),Prs(Good),ones(size(Good)),xi,yi,zi,'@nansum'),[2,1,3]);
    %relative pressure
  [xi,yi,zi] = meshgrid(Settings.Grid.Years,Settings.Grid.Months,Settings.Grid.dTP);
  Good = find(~isnan(Year+Month+dTP));
  Results.Time.dTP = Results.Time.dTP+ permute(bin2matN(3,Year(Good),Month(Good),dTP(Good),ones(size(Good)),xi,yi,zi,'@nansum'),[2,1,3]);

  clear xi yi zi Good Year Month Lon Lat TP Prs Data dTP

  
  if mod(iDay,100) == 0; save(Settings.OutFile,'Settings','Results'); end  
end; clear iDay

save(Settings.OutFile,'Settings','Results')