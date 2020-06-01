clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%generate maps of IAGOS properties (including GWs) from prepare data
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/JUN/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%general
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%variables to process
Settings.Vars = {'Prs'};%'STT_A','STT_k','T','U'};

%time period to analyse
Settings.Years = 2000:1:2010;
Settings.Days  = 1:1:30; 

%final grid size to output the results on
Settings.Grid.Lon = -180:.5:180;
Settings.Grid.Lat = -90:.5:90;

%file handling
Settings.DataDir = [LocalDataDir,'/corwin/IAGOS_annual/'];
Settings.OutFile = 'mapdata.mat';

%statistics to compute
Settings.Stats = {'mean','stdev','median','gini'};

%height bands
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pressure/height bands. Can overlap safely. Order of numbers doesn't matter.
%'t' for tropopause-relative, 'a' for absolute.
Settings.PrsBands = {{'t',-100,-25}, ...
                     {'t', -25, 25}, ...
                     {'t',  25,100}, ...
                     {'a',1000,  0}};  

%clustering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%use grid or cluster analysis?
Settings.Method = 'h'; %g for geographic map, h for hierarchical clusters

%geographic param settings only used if in geographic mode
Settings.ClusterParams.G.Lon = -180:3:180;
Settings.ClusterParams.G.Lat = -90:3:90;

%hierarchical param settings only used if in hierarchical mode
Settings.ClusterParams.H.NClusters   = 10000;
Settings.ClusterParams.H.MinPoints   = 100;
Settings.ClusterParams.H.MergeArea   = [0.25,0.25]; %points will be rounded off to the nearest this (lon/lat) and then duplicates removed before defining clusters
Settings.ClusterParams.H.MaxIter     = 1000;
Settings.ClusterParams.H.NReplicates = 10;
Settings.ClusterParams.H.MaxDist     = 500; %km from cluster centre permitted

%bootstrapping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number of samples per iteration
Settings.BS.NSamples = 1000;

%number of iterations
Settings.BS.NStraps = 100;

%number of rows to pull out of the RNG each time 
%(no effect on final results, just helps with runtime)
Settings.BS.NPerPass  = 10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create results arrays
X = NaN(numel(Settings.Grid.Lon),numel(Settings.Grid.Lat), ...
        numel(Settings.PrsBands),numel(Settings.Stats));
Results = struct();
for iVar=1:1:numel(Settings.Vars)
  Results.(Settings.Vars{iVar}) = X;
end;
Results.N   = X(:,:,:,1);
Results.Cid = X(:,:,:,1);
clear X iVar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

textprogressbar('Loading data ')
for iYear = 1:1:numel(Settings.Years);
  
  %load data for the year
  YearFile = [Settings.DataDir,'/merged_',num2str(Settings.Years(iYear)),'.mat'];
  if ~exist(YearFile,'file'); continue; end
  YearData = load(YearFile); YearData = YearData.Results;
  
  %select just the data we want
  Select.InPeriod = find(ismember(floor(date2doy(YearData.Time(:))),Settings.Days)); %on the right days
  Select.NonZero  = find((YearData.Lon + 1000.*YearData.Lat) ~= 0); %lat and lon values are not both zero
  Select.NonNaN   = find(~isnan(YearData.Lon + YearData.Lat)); %lat and lon values are not nans
  Select.GWs      = find(~isnan(YearData.STT_A)); %gravity waves were measured
  
  %combine the selections into one call
  Selected = Select.InPeriod;
  Selected = intersect(Selected,Select.NonZero);  
  Selected = intersect(Selected,Select.NonNaN); 
  Selected = intersect(Selected,Select.GWs);
  
  %pull out the indices we want
  YearData = reduce_struct(YearData,Selected);

  %store the data
  if ~exist('Store','var'); Store = YearData;
  else;                     Store = cat_struct(Store,YearData,1);
  end
  textprogressbar(iYear./numel(Settings.Years).*100)
end; 
textprogressbar('!')
clear iYear YearFile YearData Selected Select

disp('Data loaded')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% identify pressure bands, and loop over them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Bands = zeros(numel(Settings.PrsBands),numel(Store.Lon));
for iBand=1:1:numel(Settings.PrsBands);
  
  %get info
  Band = Settings.PrsBands{iBand};
  if strcmp(Band{1},'a'); Baseline = zeros(size(Store.Lon));
  else;                   Baseline = Store.TropPres;
  end
  
  %work out which bins fall in this band
  InThisBand = find((Store.Prs - Baseline) >= min(Band{2},Band{3}) ...
                  & (Store.Prs - Baseline) <= max(Band{2},Band{3}));
  
  Bands(iBand,InThisBand) = 1;
end; 
clear iBand Band InThisBand Baseline

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loop over bands
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iBand = 1:1:numel(Settings.PrsBands)
  
  disp(['Processing pressure band ',num2str(iBand),' of ',num2str(numel(Settings.PrsBands))])
  
  %select data in band
  Data = reduce_struct(Store,find(Bands(iBand,:) == 1));
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% assign data to clusters. These can be geographic or hierarchical.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  switch Settings.Method
    case 'g'; 
      %geographic assignment
      [IDs,Map] = assign_clusters_geographical(Data,Settings.Grid,Settings.ClusterParams.G); 
    case 'h'; 
      %hierarchical cluster assignment
      [IDs,Map] = assign_clusters_hierarchical(Data,Settings.Grid,Settings.ClusterParams.H); 
    otherwise; disp('Cluster assignment method not specified'); stop
  end

  %get a list of used clusters
  Clusters = unique(IDs);
  
  %only bother with filled clusters (this only applies to geographic
  %clustering - the hierarchical ones are filled by definition)
  if strcmp(Settings.Method,'g')
    Clusters = Clusters(~isnan(Clusters));
    
    %also drop them from the map
    a = unique(Map);
    a(Clusters) = [];
    for iX=1:1:numel(a); Map(Map == a(iX)) = NaN; end
    clear a iX
  end
  
  disp('--> Clusters assigned to map')  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% generate random numbers for bootstrapping
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %produce random-sampling indices (will use the same distribution for
  %each variable as this calculation is very slow)
  
  %storing all the random numbers we generate in one array is very slow
  %So, generate them in pieces and glue them together later. Go at NPerPass 
  %full samples per pass. NPerPass=10 worked well in testing.
  
  %yes, this seems dumb, but for whatever reasons (memory management??) it
  %reduces the estimated time on my test dataset from ~33 hrs to ~400s
  %I'll take some demented logic for that kind of saving
  
  textprogressbar('--> Generating randoms ')
  NPerClust = hist(IDs,Clusters);
  Randoms   = NaN(numel(Clusters),Settings.BS.NStraps,Settings.BS.NSamples,'single'); %single for memory reasons
  Piece     = NaN(numel(Clusters),Settings.BS.NPerPass,Settings.BS.NSamples,'single');
  NPieces   = ceil(Settings.BS.NStraps./Settings.BS.NPerPass);
  
  for iPiece = 1:1:NPieces;
    
    for iClust = 1:1:numel(Clusters);
      %generate randoms
      r = randi(NPerClust(iClust),10,Settings.BS.NSamples);
      
      %convert from number-within-cluster to number-within-dataset, and store
      a = find(IDs == Clusters(iClust));
      Piece(iClust,:,:) = a(r);
    end; clear iClust
    
    %store
    Randoms(:,((iPiece-1)*10)+1:iPiece*10,:) = Piece;
    
    if mod(iPiece,5) == 1; textprogressbar(iPiece./NPieces.*100); end
  end
  textprogressbar(100);textprogressbar('!')
    
  %also store the number of points associated with each cluster
  NMap = Map; NMap(:) = 0;
  for iClust = 1:1:numel(Clusters)
    NMap(Map == Clusters(iClust)) = NPerClust(iClust);
  end
  
  clear NPerClust Piece NPieces iPiece iClust r a
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% bootstrap within each cluster
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  textprogressbar('--> Bootstrapping ')
  ClusterStats = NaN(numel(Settings.Vars),numel(Clusters),numel(Settings.Stats));
  for iVar=1:1:numel(Settings.Vars);
    %get data
    V = single(Data.(Settings.Vars{iVar}));
          
    %some variable will be NaN, as not every flight measures every variable
    %replace these with points randomly sampled **from the same cluster**
    if sum(isnan(V)) > 0;
      for iCluster = 1:1:numel(Clusters);
        r = Randoms(iCluster,:,:);
        Bad  = find( isnan(V(r)));
        Good = find(~isnan(V(r)));
        if numel(Good) > 0 & numel(Bad) > 0;
          Replacement = randsample(Good,numel(Bad),'true');
          r(Bad) = Replacement;
          Randoms(iCluster,:,:) = r;
        end
      end
    end
    clear iCluster Bad Good r
    
    %compute the statistics for each cluster and bootstrap sample
    for iStat= 1:1:numel(Settings.Stats)
      switch Settings.Stats{iStat}
        case 'mean';
          Straps= mean(V(Randoms),3);
        case 'median';
          Straps = median(V(Randoms),3);
        case 'stdev';
          Straps = std(V(Randoms),[],3);
        case 'gini'
          %only defined for positive data. It's also really slow. So only use for GW amplitudes.
          %also, it only works in 2d at most. so we need to reshape it,
          %do the operation, then reshape back
          if strfind(Settings.Vars{iVar},'ST') & strfind(Settings.Vars{iVar},'_A')
            sz = size(Randoms);
            Straps= reshape(ginicoeff(reshape(V(Randoms),sz(1)*sz(2),sz(3)),2),sz(1),sz(2));
            clear sz
          end
        otherwise
          if isnumeric(Settings.Stats{iStat});
            Straps = prctile(V(Randoms),Settings.Stats{iStat},3);
          else
            disp('Statistic not specified')
          end
      end
      
      %take the median of the bootstrapped statistic
      ClusterStats(iVar,:,iStat) = nanmedian(Straps,2);
      
    end
    textprogressbar(iVar./numel(Settings.Vars).*100)    
  end
  textprogressbar('!') 
  clear iStat iVar Randoms Straps V 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% assign the results to the map
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  %put the statistics onto the map
  for iVar = 1:1:numel(Settings.Vars)
    R = Results.(Settings.Vars{iVar});
    for iStat=1:1:numel(Settings.Stats)
      Ri = R(:,:,iBand,iStat);
      for iCluster = 1:1:numel(Clusters)
        idx = find(Map' == Clusters(iCluster));
        Ri(idx) = ClusterStats(iVar,iCluster,iStat);
      end
      R(:,:,iBand,iStat) = Ri;
    end
    Results.(Settings.Vars{iVar}) = R;
  end
  
  clear iVar R iStat Ri iCluster ClusterStats Clusters idx

  %also retain the cluster IDs and the number of points per cluster
  Results.N(  :,:,iBand) = NMap'; clear NMap
  Results.Cid(:,:,iBand) = Map'; clear Map

  disp('--> Results mapped')
  
  %done!
  clear Data
end;
clear Store Bands iBand

%done!
save(Settings.OutFile,'Results','Settings')
 
 
 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% geographical cluster assignment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [IDs,Map] = assign_clusters_geographical(Data,Grid,Params)

%produce a map assigning each FINAL OUTPUT point to a cluster
[xi,yi] = meshgrid(Params.Lon,Params.Lat);
zi = xi; zi(1:1:end) = 1:1:numel(zi);
[xj,yj] = meshgrid(Grid.Lon,Grid.Lat);
Map = interp2(xi,yi,zi,xj,yj,'nearest');
clear xi yi xj yj 

%now assign all measurements to the cluster they lie within
idx_x = discretize(Data.Lon,Params.Lon);
idx_y = discretize(Data.Lat,Params.Lat);
idx = sub2ind(size(zi),idx_y,idx_x);
IDs = zi(idx);

disp('--> Data assigned to lat/lon boxes')

%done!

return; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% hierarchical cluster assignment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [IDs,Map] = assign_clusters_hierarchical(Data,Grid,Params)

%1. generate clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ok. First we need to assign each point to an initial cluster
%at the full volume of the data, this is intractable
%so, reduce the data volume by saying that any points in the same small
%area are just one point

Lon = round(Data.Lon ./ Params.MergeArea(1));
Lat = round(Data.Lat ./ Params.MergeArea(2));
[~,idx] = unique(Lat + 1000.*Lon);
clear Lat Lon

%now, generate the initial set of clusters
[~,C] = kmeans([Data.Lat(idx),Data.Lon(idx)], ...
               Params.NClusters, ...
               'Distance','sqeuclidean',...
               'MaxIter',Params.MaxIter, ...
               'Replicates',Params.NReplicates, ...
               'Display','off');

%now, re-assign all data points to a cluster, using the full lat/lon resolution.
%This is iterative, as clusters which are too underpopulated need to be removed 
%and their members reassigned

Bad = -999; %just for the first loop
while numel(Bad) > 0
  
  %find where all the un-clustered points go
  [dx,idx] = pdist2(C,[Data.Lat,Data.Lon],'euclidean','Smallest',1);
  
  %remove points too far away to be fairly included in a cluster
  idx(dx > Params.MaxDist) = [];
  
  %store the cluster ids against the data
  IDs = Data.Lon .* NaN;
  for iCluster = 1:1:size(C,1);
    IDs(idx == iCluster) = iCluster;
  end; clear iCluster
  
  %check each cluster has enough points, drop those that don't
  NPerClust = hist(IDs(:),1:1:size(C,1));
  Bad = find(NPerClust < Params.MinPoints);
  
  if numel(Bad) > 0;
    C(Bad,:) = [];
  end
end
clear NPerClust Bad dx idx

disp('--> Data assigned to hierarchical clusters')

%2. assign all geographic points to their nearest cluster centre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create results arrays
CC = NaN( numel(Grid.Lon),numel(Grid.Lat));
CD = ones(numel(Grid.Lon),numel(Grid.Lat)).*1e7; %1e7km is a fill value that should be larger than any real distance

%create lat/lon grid
[xi,yi] = ndgrid(Grid.Lon,Grid.Lat);

%now, find the closest cluster centre to each point
for iCluster = 1:1:size(C,1)

  %find distance of each point from this cluster centre
  A = [yi(:),xi(:)];
  B = A; B(:,1) = C(iCluster,1); B(:,2) = C(iCluster,2);
  dx = reshape(nph_haversine(A,B),size(xi));
  clear A B
  
  %discard any points too far from a cluster centre
  dx(dx > Params.MaxDist) = 1e7;
  
  if iCluster == 1;
    %this is the first cluster - it is the closest to every point
    CD = dx;
    CC(dx < 1e6) = 1;
  else
    %check if the cluster is closer
    Delta = CD - dx;
    Closer = find(Delta > 0);
    
    %update the distance map
    CD(Closer) = dx(Closer);
    
    %and the cluster map
    CC(Closer) = iCluster;

  end
  
end

Map = CC';

%3. final reassignment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%since we assigned each cluster to a gridpoint, it is possible (and, at high 
%cluster numbers, likely) that some clusters did not get put on the map. If
%so, find them, remove them from the ID list, and reassign the data to the 
%closest "true" cluster that actually made it onto the map

ClusterList    = unique(IDs);
ClustersMapped = unique(Map); ClustersMapped(isnan(ClustersMapped)) = [];

InBoth    = find( ismember(ClusterList,ClustersMapped));
NotInBoth = find(~ismember(ClusterList,ClustersMapped)); 

C2 = C(InBoth,:);

for iCluster = 1:1:numel(NotInBoth)
  
  InThisCluster = find(IDs == NotInBoth(iCluster));
  [~,idx] = pdist2(C2, ...
                  [Data.Lat(InThisCluster),Data.Lon(InThisCluster)], ...
                  'euclidean','Smallest',1);
  IDs(InThisCluster) = InBoth(idx);
  
end

return; end
