clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate ERA5 maps  of basic vars which can be compared to ours
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings.TimeScale = datenum(1994,8,1):1:datenum(2019,12,31);
Settings.TimeScale = datenum(2010,1,1):1:datenum(2019,12,31);
Settings.PrsRange  = [200,250];
Settings.LatScale  = 0:1.5:90-1.5;
Settings.LonScale  = -180:1.5:180-1.5;

Settings.OutFile = 'era5_iagosmap.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Results.U = NaN(numel(Settings.TimeScale),...
                numel(Settings.LonScale),numel(Settings.LatScale));
              
Results.V = Results.U;
Results.T = Results.U;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load and store data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

textprogressbar('Loading ERA5 data ')
for iDay=1:1:numel(Settings.TimeScale);
  
  %get data
  FilePath = era5_path(Settings.TimeScale(iDay));
  if ~exist(FilePath,'file'); continue; end
  Data = rCDF(FilePath);
  
  %loop over vars
  for iVar=3%1:1:3;
    
    switch iVar
      case 1; Var = Data.u;
      case 2; Var = Data.v;
      case 3; Var = Data.t;
    end
    stop
    %take daily mean
    Var = nanmean(Var,4);

    %find height range and average
    Prs = ecmwf_prs_v2([],137); idx_prs = 74;%inrange(Prs,Settings.PrsRange); clear Prs
    Var = nanmean(Var(:,:,idx_prs),3);
    clear idx_prs
    
    %find lat/lon range
    Var = Var(:,inrange(Data.latitude,[min(Settings.LatScale),max(Settings.LatScale)]));
    Var = Var(inrange(Data.longitude,[min(Settings.LonScale),max(Settings.LonScale)]),:);

    %store
    switch iVar
      case 1; Results.U(iDay,:,:) = Var;
      case 2; Results.V(iDay,:,:) = Var;
      case 3; Results.T(iDay,:,:) = Var;
    end
    
  end
  
  
  
  textprogressbar(iDay./numel(Settings.TimeScale).*100)
end
textprogressbar('!')

save(Settings.OutFile,'Settings','Results')