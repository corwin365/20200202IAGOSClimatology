clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot some example timeseries to see what is weird about them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Settings.DataDir   = [LocalDataDir,'/corwin/IAGOS_annual'];
Settings.Period    = [6,8]; %months-of-year range
Settings.MinAmp    = 1.2; %K amplitude peak in the time series
Settings.LatRange  = [60,70];
Settings.LonRange  = [-100,-60];
Settings.MinPoints = 100; %minimum number of points in the box

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load data meeting the criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for Year=1994:1:2019;
  disp(['Loading ',num2str(Year)])
  
  %load year
  ThisYear = load([Settings.DataDir,'/merged_',num2str(Year),'_v500km.mat']);
  
  %apply edgemask
  ThisYear.Results.STT_A(ThisYear.Results.STT_EdgeMask == 1) = NaN;
  
  %find months
  [~,mm,~] = datevec(ThisYear.Results.Time);
  
  %find all series which enter the box
  Flag = zeros(size(ThisYear.Results.Lat)); Flag(isnan(ThisYear.Results.Lat)) = NaN;
  
  InBox = intersect(inrange(ThisYear.Results.Lat,Settings.LatRange),...
                    inrange(ThisYear.Results.Lon,Settings.LonRange));
  InBox = intersect(InBox,find(mm >= min(Settings.Period) & mm <= min(Settings.Period)));
  Flag(InBox) = 1;
  clear InBox
  
  NFlagged = nansum(Flag,2);
  Good = find(NFlagged > Settings.MinPoints);
  clear Flag NFlagged 
  
  %trim dataset down to just these entries
  Vars = fieldnames(ThisYear.Results);
  Data = struct();
  for iVar=1:1:numel(Vars)
    V = ThisYear.Results.(Vars{iVar});
    V = V(Good,:);
    Data.(Vars{iVar}) = V;
  end; clear V iVar Vars
  
  %store
  if ~exist('Store'); Store = Data;
  else;               Store = cat_struct(Store,Data,1);
  end
  
  clear Data Good ThisYear
  
end; clear Year


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iSeries = 1:1:size(Store.Tprime,1);
  
  %find part of series in box
  InBox = intersect(inrange(Store.Lat(iSeries,:),Settings.LatRange),...
                    inrange(Store.Lon(iSeries,:),Settings.LonRange));
  
  %get data                  
  x  = 1:1:numel(Store.Tprime(iSeries,:));
  y  = Store.Tprime(iSeries,:);
  y2 = Store.STT_A( iSeries,:);
  y3 = Store.T(     iSeries,:);
  y4 = Store.U(     iSeries,:);
                  
                  
  %check peak amplitude in range
  if max(y2(InBox)) < Settings.MinAmp; continue; end
  
  clf
  subplot(4,1,1)
  plot(x,y)
  hold on
  plot(x(InBox),y(InBox));                 
  
 
  subplot(4,1,2)
  plot(x,y2)
  hold on
  plot(x(InBox),y2(InBox));   
  
 
  subplot(4,1,3)
  plot(x,y3)
  hold on
  plot(x(InBox),y3(InBox));   
    
  subplot(4,1,4)
  plot(x,y4)
  hold on
  plot(x(InBox),y4(InBox));     
  
  pause 
end