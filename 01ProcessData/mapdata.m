clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate some summary statistics for the iagos data, specifically:
%
%1. total number of cruises
%2. time distribution
%3. space distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%file handling
Settings.InDir = [LocalDataDir,'/corwin/IAGOS_st'];
Settings.OutFile = 'mapdata_200_225.mat';

%space and time gridding
Settings.Grid.t   = datenum(1994,1,1):1:datenum(1994,12,31);
Settings.Grid.Lon = -180:1:180;
Settings.Grid.Lat =  -90:1: 90;

%what pressure range are we using?
Settings.PrsRange = [200,225];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get list of files
Files = wildcardsearch(Settings.InDir,'*.mat');

%create grid for bins
Grid = zeros(numel(Settings.Grid.t), ...
             numel(Settings.Grid.Lon), ...
             numel(Settings.Grid.Lat));

[xi,yi] = meshgrid(Settings.Grid.Lon,Settings.Grid.Lat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% go!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

textprogressbar('Gridding ')
for iFile=1:1:numel(Files)
  if mod(iFile,100); textprogressbar(iFile./numel(Files).*100); end
  
  %load file
  Data = load(Files{iFile});
  
  %see if it's empty
  if nansum(Data.Results.Lat) == 0; continue; end
  
  %find the date
  FileName = Files{iFile};
  A = strfind(FileName,'ST_');
  B = strfind(FileName,'_v2.mat');
  Date = str2num(FileName(A+3:B-1));
  tidx = closest(Settings.Grid.t,Date);
  clear FileName A B Date
  
  
  %grid the lat and lon
  Lat = Data.Results.Lat(:);
  Lon = Data.Results.Lon(:);
  
  %cut down to pressire range of interest
  Prs = Data.Results.Prs(:)./100;
  InPrsRange = find(Prs > min(Settings.PrsRange) ...
                  & Prs < max(Settings.PrsRange));
  Lat = Data.Results.Lat(InPrsRange);
  Lon = Data.Results.Lon(InPrsRange);
  clear Prs InPrsRange
  
  %remove remaining empty points
  Bad = find(Lat+Lon == 0);
  Lat(Bad) = NaN;
  Lon(Bad) = NaN;
  clear Bad
  
  Good = find(~isnan(Lat+Lon));
  Lat = Lat(Good); Lon = Lon(Good);

  Grid(tidx,:,:) = bin2mat(Lon,Lat,ones(size(Lon)),xi,yi,'@nansum')';

  clear Lat Lon Good tidx
end
textprogressbar(100);textprogressbar('!')

save(Settings.OutFile,'Grid','Settings')