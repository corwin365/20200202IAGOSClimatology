clearvars


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate zonal mean plot data for IAGOS GW study
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/JUN/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Settings.DataDir     = [LocalDataDir,'/corwin/IAGOS_annual'];
Settings.TimeScale   = datenum(1994,8,1):1:datenum(2019,12,31);
Settings.AbsPrsScale = 170:10:350;
Settings.RelPrsScale = -425:5:300;
Settings.LatScale    = -40:4:80;
Settings.Vars        = {'T','U','STT_A','STT_k'};%,'STU_A','STU_k','STV_A','STV_k'};
Settings.OutFile     = 'zm_longrange.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create results arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%absolute gridded results
Results.Abs = NaN(numel(Settings.Vars),numel(Settings.TimeScale), ...
                  numel(Settings.LatScale),numel(Settings.AbsPrsScale));

%tropopause-relative gridded results                
Results.Rel = NaN(numel(Settings.Vars),numel(Settings.TimeScale), ...
                  numel(Settings.LatScale),numel(Settings.RelPrsScale));   
                
%number of points                
Results.AbsN = Results.Abs;
Results.RelN = Results.Rel;

%zonal-mean tropopause height
Results.TP  = NaN(numel(Settings.TimeScale),numel(Settings.LatScale));             
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                

Store.Name = '';
textprogressbar('Gridding data ')
for iDay=1:1:numel(Settings.TimeScale)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %get data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  [yy,~,~] = datevec(Settings.TimeScale(iDay));
  File = [Settings.DataDir,'/merged_',num2str(yy),'.mat'];
  if ~strcmp(File,Store.Name)
    Store = load(File); Store = Store.Results;
    Store.Name = File;
  end
  clear yy File
  
  OnThisDay = find(floor(Store.Time) == Settings.TimeScale(iDay));
  if numel(OnThisDay) == 0; clear OnThisDay;continue; end
  
  ThisDay = reduce_struct(Store,OnThisDay,'Name');
  clear OnThisDay
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% grid data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %get metadata needed
  Lat = ThisDay.Lat;
  AbsPrs = ThisDay.Prs;
  RelPrs = ThisDay.Prs - ThisDay.TropPres;
  
  %do binning
  for iVar=1:1:numel(Settings.Vars)

    [xi,yi] = meshgrid(Settings.LatScale,Settings.AbsPrsScale);
    Results.Abs( iVar,iDay,:,:) = bin2mat(Lat,AbsPrs,       ThisDay.(Settings.Vars{iVar}), xi,yi,'@nanmean')';
    Results.AbsN(iVar,iDay,:,:) = bin2mat(Lat,AbsPrs,~isnan(ThisDay.(Settings.Vars{iVar})),xi,yi,'@nansum')';

    [xi,yi] = meshgrid(Settings.LatScale,Settings.RelPrsScale);
    Results.Rel( iVar,iDay,:,:) = bin2mat(Lat,RelPrs,       ThisDay.(Settings.Vars{iVar}), xi,yi,'@nanmean')';
    Results.RelN(iVar,iDay,:,:) = bin2mat(Lat,RelPrs,~isnan(ThisDay.(Settings.Vars{iVar})),xi,yi,'@nansum')';

    
   if iVar == 1;
     Results.TP(iDay,:) = bin2matN(1,Lat,ThisDay.TropPres,Settings.LatScale,'@nanmean');
   end
    
  end
  clear iVar Lat AbsPrs RelPrs xi yi ThisDay
  
  if mod(iDay,500) == 0; save(Settings.OutFile,'Settings','Results','-v7.3'); end
  textprogressbar(iDay./numel(Settings.TimeScale).*100)
  
  
end
textprogressbar(100);textprogressbar('!')
clear Store iDay

save(Settings.OutFile,'Settings','Results','-v7.3')
