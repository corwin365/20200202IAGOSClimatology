clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%add tropopause height from ERA5 to each cruise from IAGOS
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/05/19
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TPSettings.DataDir.Trop  = '.';
TPSettings.DataDir.IAGOS =  [LocalDataDir,'/corwin/IAGOS_st/'];
TPSettings.TimeScale  = datenum(2012,1,1):1:datenum(2020,12,31);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%don't reload ECMWF file every year
TropData = struct();
TropData.File = '';


dd = date2doy(TPSettings.TimeScale);
[~,idx] = sort(dd,'asc');
TPSettings.TimeScale = TPSettings.TimeScale(idx);
clear idx dd


for iDay=1:1:numel(TPSettings.TimeScale)
  disp(datestr(TPSettings.TimeScale(iDay)))
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% load the tropopause data, if needed
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  [yy,~,~] = datevec(TPSettings.TimeScale(iDay));
  E5File = [TPSettings.DataDir.Trop,'/era5_tropopause_',num2str(yy),'.mat'];
  if ~strcmp(TropData.File,E5File);
    
    %load the tropopause data
    T = load(E5File);
    
    %reshape to produce a continuous timeseries for each gridpoint
    sz = size(T.Results.Tropopause);
    Tr = reshape(T.Results.Tropopause,sz(1),sz(2),sz(3)*sz(4));
    
    %create a continuous time scale
    [t,tt] = meshgrid(T.Results.t,T.Results.h);
    t = t+tt; clear tt; t = t(:);
    
    %ascending latitudes
    Lat = T.Results.Lat(end:-1:1);
    Tr = Tr(:,end:-1:1,:);
    
    %create a global interpolant for this year
    TropData.I = griddedInterpolant({T.Results.Lon,Lat,t},Tr);
    
    %store the name of the tile so we don't do this again
    TropData.File  = E5File;
    
    %tidy up
    clear T sz Tr t Lat
    
  end
  clear yy E5File

  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% load the daily STed IAGOS file and
  % interpolate the tropopause onto all the tracks
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %load data
  DayFile = wildcardsearch(TPSettings.DataDir.IAGOS,['*',num2str(TPSettings.TimeScale(iDay)),'*']);
  if numel(DayFile) == 0; clear DayFile; continue; end
  Data = load(DayFile{1});
  
  %interpolate data
  Data.Results.TropPres = single(TropData.I(Data.Results.Lon,Data.Results.Lat,Data.Results.Time));
  
  
  %save
  Settings = Data.Settings;
  Results  = Data.Results;
  save(DayFile{1},'Settings','Results')
  clear DayFile Settings Results Data
  
  
  
  
  
  
  
end