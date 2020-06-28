clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate individual time series for each identified peak
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time and pressure range
Settings.TimeScale = datenum(1994,1,1):31:datenum(2019,12,31);
Settings.PrsRange = [10000,0];

%how much to smooth in time?
Settings.SmoothSize = 3;

%which variables to use
Settings.Vars = {'STT_A','STT_k','U','T'};

%which cluster map do we want to base our analysis on?
Settings.Cluster.Dir = '../02Maps/out/';
Settings.Cluster.Mode = 'h';
Settings.Cluster.Layer = 1;
Settings.Cluster.Suffix = 'sgolay900';
Settings.Cluster.Quarter = 'DJF'; 

%regression settings
Settings.Regression.Lags = -11:1:11;
Settings.Regression.Indices = {'QBO','ENSO','HadCRUT','NAM','TSI'};

% % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % %% data import
% % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % 
% % % % % % %first, load the cluster map 
% % % % % % Data = load([Settings.Cluster.Dir,'/',Settings.Cluster.Mode, ...
% % % % % %              '_',Settings.Cluster.Quarter,'_b',num2str(Settings.Cluster.Layer), ...
% % % % % %              '_',Settings.Cluster.Suffix,'.mat']);
% % % % % % Map.Cid = Data.Results.Cid;
% % % % % % Map.Lon = Data.Settings.Grid.Lon;
% % % % % % Map.Lat = Data.Settings.Grid.Lat;
% % % % % % clear Data
% % % % % % 
% % % % % % % % %temporary for testing - trim to newfoundland-ish
% % % % % % % % lon = inrange(Map.Lon,[-70,-50]);
% % % % % % % % lat = inrange(Map.Lat,[40,50]);
% % % % % % % % Map.Lon = Map.Lon(lon); Map.Lat = Map.Lat(lat); Map.Cid = Map.Cid(lon,:); Map.Cid = Map.Cid(:,lat);
% % % % % % % % clear lon lat
% % % % % % 
% % % % % %   
% % % % % % 
% % % % % % %second, load the aircraft data
% % % % % % [yy1,~,~] =datevec(min(Settings.TimeScale));
% % % % % % [yy2,~,~] =datevec(max(Settings.TimeScale));
% % % % % % clear Store; %just to be certain, as needs to be gone for internal logic
% % % % % % for Year=yy1:1:yy2;
% % % % % %   disp(['Loading ',num2str(Year)])
% % % % % %   
% % % % % %   %load storage file
% % % % % %   Data = load([LocalDataDir,'/corwin/IAGOS_annual/merged_',num2str(Year),'_sgolay900.mat']);
% % % % % %  
% % % % % %   %drop bad points
% % % % % %   Good = find((1000.*Data.Results.Lon+Data.Results.Lat) ~= 0 ...
% % % % % %             & 1./Data.Results.STT_k > 25             ...
% % % % % %             & ~isnan(Data.Results.Lon + Data.Results.Lat + Data.Results.STT_A));
% % % % % %   Data = reduce_struct(Data.Results,Good);
% % % % % %   
% % % % % %   %pull out what we need, and store
% % % % % %   if ~exist('Store'); Store = Data; else  Store = cat_struct(Store,Data,1); end
% % % % % % 
% % % % % %   
% % % % % % end; clear yy1 yy2 Year Good Data
% % % % % % 
% % % % % % 
% % % % % % 
% % % % % % 
% % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % %% produce an individual monthly time series for every cluster
% % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % 
% % % % % % [FirstYear,FirstMonth,~] = datevec(min(Settings.TimeScale));
% % % % % % [ LastYear, LastMonth,~] = datevec(max(Settings.TimeScale));
% % % % % % MonthScale = [];
% % % % % % for iY=FirstYear:1:LastYear
% % % % % %   if iY == FirstYear;   m1 = FirstMonth; else FirstMonth = 1; end
% % % % % %   if iY == LastYear;  mend =  LastMonth; else mend = 12; end
% % % % % %   for iM=m1:1:mend; MonthScale(end+1) = iM + 12.*(iY-FirstYear); end
% % % % % % end; clear iY m1 mend iM FirstMonth LastYear LastMonth
% % % % % % 
% % % % % % 
% % % % % % ClusterSeries = NaN(numel(Settings.Vars),numel(MonthScale),nanmax(Map.Cid(:)));
% % % % % % 
% % % % % % textprogressbar('Producing monthly cluster data ')
% % % % % % for iMonth=1:1:numel(MonthScale)
% % % % % % 
% % % % % %   %find points in month
% % % % % %   ThisMonth = inrange(Store.Time, ...
% % % % % %                       [datenum(FirstYear,MonthScale(iMonth),1), ...
% % % % % %                        datenum(FirstYear,MonthScale(iMonth)+1,-1)]);
% % % % % %   if numel(ThisMonth) == 0; continue; end
% % % % % %   
% % % % % %   %take the data and identify which cluster it falls into
% % % % % %   Store.Cid = NaN(numel(ThisMonth),1);
% % % % % %   for iPoint = 1:1:numel(ThisMonth)
% % % % % %     idx_lat           = closest(Map.Lat,Store.Lat(ThisMonth(iPoint)));
% % % % % %     idx_lon           = closest(Map.Lon,Store.Lon(ThisMonth(iPoint)));
% % % % % %     Store.Cid(iPoint) = Map.Cid(idx_lon,idx_lat);
% % % % % %   end; clear iPoint idx_lon idx_lat
% % % % % %   
% % % % % %   %hence, find the monthly mean for each cluster
% % % % % %   for iVar=1:1:numel(Settings.Vars)
% % % % % %     Var = Store.(Settings.Vars{iVar});
% % % % % %     Var = Var(ThisMonth);
% % % % % %     for iCluster=1:1:nanmax(Map.Cid(:))
% % % % % %       ClusterSeries(iVar,iMonth,iCluster) = nanmedian(Var(Store.Cid == iCluster));
% % % % % %     end
% % % % % %   end
% % % % % %   clear iVar Var iCluster
% % % % % %   textprogressbar(iMonth./numel(MonthScale).*100)
% % % % % % end
% % % % % % textprogressbar(100); textprogressbar('!')
% % % % % % clear iMonth ThisMonth Store
% % % % % % TimeScale = datenum(FirstYear,MonthScale,15); clear FirstYear MonthScale
% % % % % % 
% % % % % % 
% % % % % % 
% % % % % % 
% % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % %% load indices
% % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % 
% % % % % % %load the indices, and put on the same monthly scale
% % % % % % Indices = NaN(numel(Settings.Regression.Indices),numel(TimeScale));
% % % % % % for iIndex=1:1:numel(Settings.Regression.Indices)
% % % % % %   
% % % % % %   switch Settings.Regression.Indices{iIndex}
% % % % % %     case 'QBO'
% % % % % %       QBO = load([LocalDataDir,'/Miscellany/QBO.mat']);
% % % % % %       QBO.Time = floor(QBO.Time); %shift from noon to midnight to make the logic easier - on a 91-day smoothing this is very minor...
% % % % % %       a = interp1(QBO.Time,QBO.QBO,TimeScale);
% % % % % %       clear QBO
% % % % % %     case 'ENSO'
% % % % % %       ENSO = load([LocalDataDir,'/Miscellany/nino34.mat']);
% % % % % %       a = interp1(ENSO.Time,ENSO.Nino34,TimeScale);
% % % % % %       clear ENSO
% % % % % %     case 'HadCRUT'
% % % % % %       HadCRUT = rCDF([LocalDataDir,'/Miscellany/HadCRUT.4.6.0.0.median.nc']);
% % % % % %       HadCRUT.MatlabTime = datenum(1850,1,HadCRUT.time);
% % % % % %       HadCRUT.NH = squeeze(nanmean(HadCRUT.temperature_anomaly(:,HadCRUT.latitude > 0,:),[1,2]));
% % % % % %       a = interp1(HadCRUT.MatlabTime,HadCRUT.NH,TimeScale);
% % % % % %       clear HadCRUT
% % % % % %     case 'NAM'
% % % % % %       NAM = load([LocalDataDir,'/Miscellany/daily_nam.mat']);
% % % % % %       a = interp1(NAM.Time,NAM.NAM,TimeScale);
% % % % % %       clear NAM
% % % % % %     case 'TSI'
% % % % % %       TSI = load([LocalDataDir,'/Miscellany/tsi.mat']);
% % % % % %       a = interp1(TSI.Time,TSI.TSI,TimeScale);
% % % % % %       clear TSI
% % % % % %   end
% % % % % %   
% % % % % %   Indices(iIndex,:) = a;
% % % % % %   clear a
% % % % % %   
% % % % % % end; clear iIndex
% % % % % % 
% % % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% regression time!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load reg.mat
Settings.SmoothSize = 3;

%smooth the indices (data is done inside the loop)
Indices = smoothn2(Indices,[1,Settings.SmoothSize]);

%create storage array
Reg.Est = NaN(numel(Settings.Vars),nanmax(Map.Cid(:)),numel(Settings.Regression.Lags),size(Indices,1));
Reg.SE  = Reg.Est;
Reg.T   = Reg.Est;
Reg.P   = Reg.Est;

textprogressbar('Regressing clusters ')
for iCluster=1:1:nanmax(Map.Cid(:));
  for iVar=1:1:numel(Settings.Vars)

    %get the data
    TimeSeries = ClusterSeries(iVar,:,iCluster);
    if nansum(TimeSeries) == 0;continue; end
    
    for iLag = 1:1:numel(Settings.Regression.Lags)
      
      %lag the data
      Lag = Settings.Regression.Lags(iLag);
      A = circshift(TimeSeries,Lag);
      
      %remove points which have gone off the end
      if     Lag > 0; A(      1:Lag) = NaN;
      elseif Lag < 0  A(end-Lag:end) = NaN;
      end
      
      %smooth
      A = smoothn2(A,Settings.SmoothSize);
      
      %do the regression
      warning off
      mdl = fitlm(Indices',A);
      warning on
      Coef = table2array(mdl.Coefficients);
      
      %and store the outputs
      Reg.Est(iVar,iCluster,iLag,:) = Coef(2:end,1);
      Reg.SE( iVar,iCluster,iLag,:) = Coef(2:end,2);
      Reg.T(  iVar,iCluster,iLag,:) = Coef(2:end,3);
      Reg.P(  iVar,iCluster,iLag,:) = Coef(2:end,4);
      
      
    end
  end
  textprogressbar(iCluster./nanmax(Map.Cid(:)).*100)
end
clear iCluster iVar TimeSeries iLag Lag A Coef mdl
textprogressbar(100); textprogressbar('!')
clear iVar iCluster iLag



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find the best coefficients, and put them on the map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,Best] = nanmax(abs(Reg.Est),[],3);

RegMaps = NaN([numel(Settings.Vars),numel(Settings.Regression.Indices),numel(Map.Cid),4]);
textprogressbar('Mapping regression outputs ')
for iVar=1:1:numel(Settings.Vars)
  for iIndex=1:1:numel(Settings.Regression.Indices)
    
    %find optimal value for each cluster
    for iCluster=1:1:nanmax(Map.Cid(:))
      Est = squeeze(Reg.Est(iVar,:,Best(iVar,iCluster,1,iIndex),iIndex));
      SE  = squeeze(Reg.SE( iVar,:,Best(iVar,iCluster,1,iIndex),iIndex));
      T   = squeeze(Reg.T(  iVar,:,Best(iVar,iCluster,1,iIndex),iIndex));
      P   = squeeze(Reg.P(  iVar,:,Best(iVar,iCluster,1,iIndex),iIndex));
    end
    
    %put on map
    for iCluster=1:1:nanmax(Map.Cid(:))
      ThisCluster = find(Map.Cid == iCluster);
      RegMaps(iVar,iIndex,ThisCluster,1) = Est(iCluster);
      RegMaps(iVar,iIndex,ThisCluster,2) = SE( iCluster);
      RegMaps(iVar,iIndex,ThisCluster,3) = T(  iCluster);
      RegMaps(iVar,iIndex,ThisCluster,4) = P(  iCluster);
    end
          
  end
textprogressbar(iVar./numel(Settings.Vars).*100)
end
textprogressbar(100); textprogressbar('!')
clear Best iVar iIndex iCluster ThisCluster
sz = size(RegMaps);
RegMaps = reshape(RegMaps,sz(1),sz(2),numel(Map.Lon),numel(Map.Lat),sz(4));
clear sz

% save('reg.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% and save to plot elsewhere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m_proj('lambert','lat',[25,80],'lon',[-130 150])
Data = squeeze(RegMaps(1,2,:,:,1)); m_pcolor(Map.Lon,Map.Lat,smoothn2(Data,[1,1].*15)'); caxis(prctile(Data(:),[10,90]));; redyellowblue16
m_coast('color','k')
m_grid