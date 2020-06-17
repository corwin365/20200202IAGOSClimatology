clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find annual cycle of UTLS wind variance in Aeolus
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/06/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%data paths
Settings.DataDir = [LocalDataDir,'/Aeolus/NC_FullQC/'];
Settings.OutFile = 'aeolus_north.mat';

%Aeolus formatting
Settings.Format  = '2B';

%data gridding
Settings.TimeScale  = datenum(2020,5,01):1:datenum(2020,06,1);
Settings.LatScale   = 40:10:90;
Settings.LonScale   = -130:20:150;
Settings.HeightLevs = 5:2.5:15;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%measurement grid
[xi,yi,zi] = meshgrid(Settings.LonScale,Settings.LatScale,Settings.HeightLevs);

%results array
Results = NaN([numel(Settings.TimeScale),size(xi)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load, grid, store
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

textprogressbar('Gridding stdevs ')
for iDay=1:1:numel(Settings.TimeScale)
  
  
  %get data for this day
  [yy,mm,dd] = datevec(Settings.TimeScale(iDay));
  OnThisDay = wildcardsearch(Settings.DataDir, ...
                             ['*AE_',Settings.Format,'_', ...
                             sprintf('%04d',yy),'-',sprintf('%02d',mm),'-', ...
                             sprintf('%02d',dd),'*']);
  clear yy mm dd
  if numel(OnThisDay) == 0; clear OnThisDay; continue; end
  
  %loop over files and store the data
  clear Store
  for iFile=1:1:numel(OnThisDay)
    
    %load data
    Data = rCDF(OnThisDay{iFile});
    
    %apply QC flags and filter controls
    Bad = find(Data.QC_Flag_Both == 0);
    Data.Zonal_wind_projection(Bad) = NaN;
    Bad = find(Data.Zonal_wind_projection < -20000);
    Data.Zonal_wind_projection(Bad) = NaN;
    
    %fix lons into 180:180
    Data.lon(Data.lon > 180) = Data.lon(Data.lon > 180) - 360;
    
    %convert alt to km
    Data.alt = Data.alt./1000;
        
    %store data
    if ~exist('Store','var'); Store = Data;
    else;                     Store = cat_struct(Store,Data,1,{'MetaData','RG'});
    end
    
  end; clear iFile Data Bad  OnThisDay
  
  %find variance of zonal wind on our output grid
  Good = find(~isnan(Store.lon + Store.lat + Store.alt + Store.Rayleigh_HLOS_wind_speed));
  zz = bin2matN(3,Store.lon(Good),Store.lat(Good),Store.alt(Good), ...
                  Store.Rayleigh_HLOS_wind_speed(Good), ...
                  xi,yi,zi,'@nanstd');
  %store
  Results(iDay,:,:,:) = zz;
  
  %tidy
  clear Good zz
  
  textprogressbar(iDay./numel(Settings.TimeScale).*100)
end; clear xi yi zi iDay
textprogressbar('!')

save(Settings.OutFile,'Settings','Results')