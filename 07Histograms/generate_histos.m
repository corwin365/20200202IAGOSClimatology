clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%generate data for histograms of wave properties in IAGOS data
%Corwin Wright, c.wright@bath.ac.uk, 23/MAY/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Settings.DataDir   = [LocalDataDir,'/corwin/IAGOS_annual'];
Settings.TimeScale = datenum(1994,1,1):1:datenum(2019,1,365);


%variables we want, and the bin sizes for them
Settings.Vars.STT_A = logspace(-3,3,1000);
Settings.Vars.STT_k = logspace(-4,-2,100);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% results storage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Results = struct();

Vars = fieldnames(Settings.Vars);
for iVar=1:1:numel(Vars)
  Results.(Vars{iVar}) = zeros(numel(Settings.Vars.(Vars{iVar})),1);
end; clear iVar Vars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% go!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AllData.Name= '' ;%avoid duplicate file loading

for iDay=1:1:numel(Settings.TimeScale)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% load data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %load file if needed
  [yy,~,~] = datevec(Settings.TimeScale(iDay));
  File = [Settings.DataDir,'/merged_',num2str(yy),'.mat']; 
  if strcmp(File,AllData.Name) ~= 1;
    if ~exist(File,'file'); clear File; continue; end
    AllData = load(File); AllData = AllData.Results;
    AllData.Name = File; clear File
    AllData.Time(AllData.Time == 0) = NaN;
    disp(yy)
  end; clear File yy
    
  
  %pull out this day
  OnThisDay = find(floor(nanmean(AllData.Time,2)) == Settings.TimeScale(iDay));
  if numel(OnThisDay) == 0; clear OnThisDay; continue; end
  Vars = fieldnames(AllData); Data = struct();
  for iVar = 1:1:numel(Vars)
    if strcmp(Vars{iVar},'Name'); continue; end
    V = AllData.(Vars{iVar});
    Data.(Vars{iVar}) = V(OnThisDay,:);
  end; clear iVar Vars V OnThisDay
  
  %take histogram of each variable we want
  Vars = fieldnames(Settings.Vars);
  for iVar=1:1:numel(Vars)
    D = Data.(Vars{iVar});
    Results.(Vars{iVar}) = Results.(Vars{iVar}) + hist(D(:),Settings.Vars.(Vars{iVar}))';
  end
  
  
end