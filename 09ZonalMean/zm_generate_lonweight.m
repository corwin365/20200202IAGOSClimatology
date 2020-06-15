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
Settings.AbsPrsScale = 165:10:350;
Settings.RelPrsScale = -425:500:300; %i.e. keep the code but run as fast as possible
Settings.LatScale    = -40:4:80;
Settings.LonScale    = -180:4:180; %used for weighting, not in final product
Settings.Vars        = {'U','STT_A'};%,'STT_k'};%'T','STU_A','STU_k','STV_A','STV_k'};
Settings.OutFile     = 'zm_final_lonweighted.mat';

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
  File = [Settings.DataDir,'/merged_',num2str(yy),'_v7.mat'];
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
  Lon = ThisDay.Lon;
  AbsPrs = ThisDay.Prs;
  RelPrs = ThisDay.Prs - ThisDay.TropPres;
  
  %do binning
  for iVar=1:1:numel(Settings.Vars)

    [xi,yi,zi] = meshgrid(Settings.LonScale,Settings.LatScale,Settings.AbsPrsScale);
    Results.Abs( iVar,iDay,:,:) = squeeze(nanmean(bin2matN(3,Lon,Lat,AbsPrs,ThisDay.(Settings.Vars{iVar}), xi,yi,zi,'@nanmean'),2));

    [yi,zi] = meshgrid(Settings.LatScale,Settings.AbsPrsScale);
    Results.AbsN(iVar,iDay,:,:) = bin2mat(Lat,AbsPrs,~isnan(ThisDay.(Settings.Vars{iVar})),yi,zi,'@nansum')';

    [xi,yi,zi] = meshgrid(Settings.LonScale,Settings.LatScale,Settings.RelPrsScale);
    Results.FRel( iVar,iDay,:,:) = squeeze(nanmean(bin2matN(3,Lon,Lat,RelPrs,ThisDay.(Settings.Vars{iVar}), xi,yi,zi,'@nanmean'),2));    
    
    [yi,zi] = meshgrid(Settings.LatScale,Settings.RelPrsScale);    
    Results.RelN(iVar,iDay,:,:) = bin2mat(Lat,RelPrs,~isnan(ThisDay.(Settings.Vars{iVar})),yi,zi,'@nansum')';

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
