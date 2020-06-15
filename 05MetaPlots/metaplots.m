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
Settings.DataDir = [LocalDataDir,'/corwin/IAGOS_annual'];
Settings.OutFile = 'metadata_all_v7.mat';


%map data
Settings.Grid.Lon = -180:5:180;
Settings.Grid.Lat = -90:5:90;
Settings.Grid.dTP = -150:5:150; %tropopause relative pressure (hPa)
Settings.Grid.Prs = 150:10:350;  %absolute pressure (hPa)

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

%make a list of months to loop over
[yyi,mmi] = meshgrid(Settings.Grid.Years,Settings.Grid.Months);
yyi = yyi(:); mmi = mmi(:);

AllData.Name = ''; %avoid reloading data unnecessarily


for iMonth= 1:1:numel(mmi)
  
  if mmi(iMonth) == 1;
    disp(['Processing ',num2str(yyi(iMonth))])
  end
  
  %get data
  File = [Settings.DataDir,'/merged_',num2str(yyi(iMonth)),'_v7.mat'];
  if strcmp(AllData.Name,File) ~= 1;
    %load the data
    if ~exist(File,'file'); clear File; continue; end
    AllData = load(File); 
    AllData.Name = File; clear File
    AllData.Results.Time(AllData.Results.Time == 0) = NaN;
    AllData.Day = floor(nanmean(AllData.Results.Time,2));
    [AllData.yy,AllData.mm,~]= datevec(AllData.Day);
  end
  
  %pull out just the month we want
  ThisMonth = find(AllData.yy == yyi(iMonth) & AllData.mm == mmi(iMonth));
  if numel(ThisMonth) == 0; clear OnThisDay;continue; end
  Data = AllData;
  Vars = fieldnames(Data.Results);
  for iVar=1:1:numel(Vars);
    V = Data.Results.(Vars{iVar});
    V = V(ThisMonth,:);
    Data.Results.(Vars{iVar}) = V;
  end; clear V iVar Vars;

  
  %pull out data we need, and tidy it up
  Good = find(~isnan(Data.Results.STT_A) & Data.Results.STT_A~= 0 & 1./Data.Results.STT_k > 80);
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

  
end; clear iDay

save(Settings.OutFile,'Settings','Results');disp('Saved');