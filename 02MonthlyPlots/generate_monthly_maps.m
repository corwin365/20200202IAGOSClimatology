clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate monthly maps of IAGOS postprocessed data, composited over all years
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/05/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%file handling
Settings.DataDir = [LocalDataDir,'/corwin/IAGOS_ST/'];
Settings.OutFile = 'IAGOS_maps_5deg.mat';

%gridding
Settings.Lon = -180:5:180;
Settings.Lat =-90:5:90;

%variables
Settings.Vars = {'A','k','Prs','Z','U','V','T'};

%outlier cutoff
Settings.CutOff = [5,95];

%years to use
Settings.Years = 1994:1:2020;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create results arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = NaN(12,numel(Settings.Lon),numel(Settings.Lat));
Results = struct();
for iVar=1:1:numel(Settings.Vars)
  Results.(Settings.Vars{iVar}) = X;
end;
clear iVar X

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do it!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%loop over months
for iMonth = 1:1:12 %hopefully this is an uncontentious number of months
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % get data for this month over all years
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  %find all days in this month over all years
  Days = [];
  for  iYear=1:1:numel(Settings.Years)
    Days = [Days,datenum(Settings.Years(iYear),iMonth,1):1:datenum(Settings.Years(iYear),iMonth+1,1)-1];
  end
  
  %create temporary data storage arrays
  Data = struct();
  Vars = [Settings.Vars,'Lat','Lon']; %geoloc needed for binning
  for iVar=1:1:numel(Vars);
    Data.(Vars{iVar}) = [];
  end
  
  %load data
  textprogressbar(['Loading data for ',datestr(Days(1),'mmmm'),' ']);
  for iDay=1:1:numel(Days)
    textprogressbar(iDay./numel(Days).*100);
    
    %find file for this day
    ThisDayFile = wildcardsearch(Settings.DataDir,['*',num2str(Days(iDay)),'*']);
    if numel(ThisDayFile) == 0; clear ThisDayFile; continue; end
    
    %load data
    ThisDayData = load(ThisDayFile{1});
    clear ThisDayFile
    
    %store what we need
    for iVar=1:1:numel(Vars);
      V = ThisDayData.Results.(Vars{iVar});
      Data.(Vars{iVar}) = cat(1,Data.(Vars{iVar}),V(:));
    end

    %and done!
    clear ThisDayData iVar V
    
  end
  textprogressbar('!')
  
  clear Vars % we DO NOT want to later confuse this with Settings.Vars, as it includes geolocation
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % map and store the data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
  textprogressbar('Binning...')
  for iVar=1:1:numel(Settings.Vars);
    %get data
    V = Data.(Settings.Vars{iVar});
    
    %remove outliers
    Bad = prctile(V,Settings.CutOff);
    V(V < min(Bad)) = NaN;
    V(V > max(Bad)) = NaN;
    
    %grid
    [xi,yi] = meshgrid(Settings.Lon,Settings.Lat);
    zz = bin2matN(2,Data.Lon,Data.Lat,V,xi,yi,'@nanmean');

    %store
    R = Results.(Settings.Vars{iVar});
    R(iMonth,:,:) = zz';
    Results.(Settings.Vars{iVar}) = R;
    
    clear V xi yi zz R Bad
  end; clear iVar  
  textprogressbar(' done')
  
  %save
  save(Settings.OutFile,'Results','Settings')
  
end; clear iMonth

%done. save
save(Settings.OutFile,'Results','Settings')