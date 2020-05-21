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
% TPSettings.TimeScale  = datenum(1994,8,1):1:datenum(2019,12,31);
TPSettings.TimeScale  = [datenum(2002,1,1):1:datenum(2005,12,31),datenum(2009,1,1):1:datenum(2011,12,31)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%don't reload ECMWF file every year
TropData = struct();
TropData.File = '';



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
  %also do some maintenance and sanity-checking
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %load data
  DayFile = wildcardsearch(TPSettings.DataDir.IAGOS,['*',num2str(TPSettings.TimeScale(iDay)),'*']);
  if numel(DayFile) == 0; clear DayFile; continue; end
  Data = load(DayFile{1});
  
  %interpolate data
  Data.Results.TropPres = single(TropData.I(Data.Results.Lon,Data.Results.Lat,Data.Results.Time));
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% also do some maintenance and sanity-checking
  %
  %all of this should semantically be in the routine 
  %that generates the data, but was identified after
  %running it for several days to make all the data, 
  %so I didn't want to run it all again with this 
  %code in the logically-correct place. It works fine 
  %here, it's just inelegant.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  Vars = fieldnames(Data.Results);
  
  %remove points where both lat and lon equal zero
  Bad = find(Data.Results.Lat + 1000.*Data.Results.Lon == 0);
  if numel(Bad) > 0;
    for iVar = 1:1:numel(Vars);
      V = Data.Results.(Vars{iVar});
      V(Bad) = NaN;
      Data.Results.(Vars{iVar}) = V;
    end; clear iVar V 
  end; clear Bad
    
    
  %drop empty columns (left in by parent routine in certain edge cases,
  %now is as good a time as any to remove)
  a = max(find(nansum(Data.Results.Lat,1) ~= 0)); %if the longest cruise in the file ends *exactly* on the equator and is going northwards, this will drop the last point. this is unlikely.
  if a < size(Data.Results.STT_A,2)
      for iVar = 1:1:numel(Vars); 
        V = Data.Results.(Vars{iVar}); 
        V = V(:,1:max(a));
        Data.Results.(Vars{iVar}) = V;
      end

  end
  clear a iVar
  
  %some time series crossing the dateline have continuous data, but with a
  %large lump of NaNs in the middle. Not sure what is causing these, and
  %the number of flights is so small it's not worth rerunning the whole
  %dataset to fix. So just find and fix them. But be careful,and do lots of
  %testing to make sure these are the droids we're looking for.
  for iCruise = 1:1:size(Data.Results.Lon,1);
    Cruise.Lon = Data.Results.Lon(iCruise,:);
    Jump = find(diff(find(~isnan(Cruise.Lon))) > 1);
    
    if numel(Jump) ~= 1; continue; end %this ony occurs once in such records
    if abs(Cruise.Lon(Jump)) < 175; continue; end %must be near dateline
    Next = min(find(~isnan(Cruise.Lon(Jump+2:end))));
    if abs(Cruise.Lon(Jump+Next+1)) < 175; continue; end %must be near dateline  
    
    %ok, stitch the time series together.
    Good = [1:Jump,Jump+Next+1:numel(Cruise.Lon)];
    Order = unique([Good,1:1:numel(Cruise.Lon)],'stable');
    for iVar = 1:1:numel(Vars);
      V = Data.Results.(Vars{iVar});
      V(iCruise,:) = V(iCruise,Order);
      Data.Results.(Vars{iVar}) = V;
    end; clear iVar V Good Order
    
    %check total length again, and truncate if needed
    a = max(find(nansum(Data.Results.Lat,1) ~= 0)); %if the longest cruise in the file ends *exactly* on the equator and is going northwards, this will drop the last point. this is unlikely.
    if a < size(Data.Results.STT_A,2)
      for iVar = 1:1:numel(Vars);
        V = Data.Results.(Vars{iVar});
        V = V(:,1:max(a));
        Data.Results.(Vars{iVar}) = V;
      end
      
    end
    clear a iVar

  end
  
  clear Vars
  
    
  %there is some inconsistency in the input data about the units of
  %pressure, which I did not account for when generating the data. Most of
  %each flight should be at height levels in the hundreds of hPa, so check
  %that the mean is of this order or magnitude and scale if not
  fac = round(log10(nanmean(round(Data.Results.Prs(Data.Results.Prs ~= 0).*100)./100)));
  if  fac == 2;
%       disp('PRESSURES FINE')
  else
%       disp(['CORRECTING PRESSURES BY FACTOR 10^',num2str(fac)])
    Data.Results.Prs = Data.Results.Prs./ 10.^(fac-2);
%     disp(['Mean now ',num2str(round(nanmean(Data.Results.Prs(Data.Results.Prs ~= 0))))])
  end
  clear fac
  

  
  %hand check any records longer than 10 000km, as these shouldn't exist,
  %and manually fix them.
  %(the steps in this section above were produced by examining stops here,
  %so this should not actually fire unless more data is acquired with more weird foibles)
  if size(Data.Results.Lat,2) > 5000; stop; end

  Settings = Data.Settings;
  Results  = Data.Results;
  save(DayFile{1},'Settings','Results')
  clear DayFile Settings Results Data
  
  
  
  
  
  
  
end