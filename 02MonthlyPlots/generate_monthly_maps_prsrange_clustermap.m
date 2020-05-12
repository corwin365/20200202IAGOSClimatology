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
Settings.DataDir = [LocalDataDir,'/corwin/IAGOS_st/'];
Settings.OutFile = 'clustertest.mat';

%gridding - this is the final map gridding, based on where the clusters
%fall, and does not represent the area averaged over.
Settings.Lon = -180:1:180;
Settings.Lat = -90:1:90;

%restrict pressure range
Settings.PrsRange = [200,250]; %hPa

%variables
Settings.Vars = {'A','k','Prs','Z','U','V','T'};

%number of clusters
Settings.NClusters = 5000;

%maximum distance of a point from cluster centre (km)
Settings.MaxDist = 500;

%years to use
Settings.Years = 1994:1:2020;%;1994:1:2000;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create results arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = NaN(12,numel(Settings.Lon),numel(Settings.Lat));
Results = struct();
for iVar=1:1:numel(Settings.Vars)
  Results.(Settings.Vars{iVar}) = X;
end;
clear iVar X

Results.Cid = Results.(Settings.Vars{1}); %cluster id
Results.N   = Results.(Settings.Vars{1}); %number of points actually used in cluster

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
    
    %discard any data outside our pressure range
    InPrsRange = find(ThisDayData.Results.Prs./100 >= min(Settings.PrsRange) ...
                    | ThisDayData.Results.Prs./100 <= max(Settings.PrsRange));
   
    %store what we need
    for iVar=1:1:numel(Vars);
      V = ThisDayData.Results.(Vars{iVar});
      Data.(Vars{iVar}) = cat(1,Data.(Vars{iVar}),flatten(V(InPrsRange)));
    end


    %and done!
    clear ThisDayData iVar V InPrsRange
    
  end
  textprogressbar('!')
  

  
  %my above algorithm leaves some zeros at the end where a longer time
  %series was added to a smaller one. these will all have [lat,lon] = [0,0]. 
  %remove these. We'll also lose any real data at exactly [0,0], but this
  %is pretty rare so it might be one or two points in 20 years of data.

  Bad = find(Data.Lon == 0 & Data.Lat == 0);
  Good = 1:1:numel(Data.Lon); Good(Bad) = [];
  for iVar=1:1:numel(Vars);
      V = Data.(Vars{iVar});
      Data.(Vars{iVar}) = V(Good);
  end
    
  %also remove NaNs, to save time later
  Good = find(~isnan(Data.Lon) & ~isnan(Data.Lat));
  for iVar=1:1:numel(Vars);
      V = Data.(Vars{iVar});
      Data.(Vars{iVar}) = V(Good);
  end
  clear iVar V Good Bad
  
  clear Vars % we DO NOT want to later confuse this with Settings.Vars, as it includes geolocation    
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % OK. we have the data lists. Now we need to make the maps
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %k-means cluster the data spatially
  %do this in lat and lon coordinates, as spatial are just as bad due to
  %the large area covered
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %we have far too many points to be tractable in computer time
  %but we know that they're very heavily clumped, since aircraft follow
  %fixed routes. So, let's treat anything in a .2x.2 degree box as
  %'saturated' and as a single point for the cluster generation process
  Lat = round(Data.Lat.*5);
  Lon = round(Data.Lon.*5);
  Combo = Lat + 1000.*Lon;
  [~,idx] = unique(Combo);
  clear Combo Lat Lon
  
  stream = RandStream('mlfg6331_64');  % Random number stream
  options = statset('UseParallel',1,'UseSubstreams',1,...
                    'Streams',stream);
  [idx,C] = kmeans([Data.Lat(idx),Data.Lon(idx)], ...
                   Settings.NClusters, ...
                   'Distance','cityblock',...
                   'MaxIter',1000, ...
                   'Replicates',10, ...
                   'Display','final', ...
                   'Options',options);


  %now, find where all the un-clustered points go
  [~,idx] = pdist2(C,[Data.Lat,Data.Lon],'euclidean','Smallest',1);
  
  %store the cluster ids against the data
  Data.Cluster = Data.Lon .* NaN;
  for iCluster = 1:1:size(C,1);
    Data.Cluster(idx == iCluster) = iCluster;
  end; clear iCluster
 
 
  clear idx
  
  %% assign all geographic points to their nearest cluster centre
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %create results arrays
  CC = NaN( numel(Settings.Lon),numel(Settings.Lat));
  CD = ones(numel(Settings.Lon),numel(Settings.Lat)).*1e6;
  
  %create lat/lon grid
  [xi,yi] = ndgrid(Settings.Lon,Settings.Lat);

  %now, find the closest cluster centre to each point
  for iCluster = 1:1:size(C,1)
    

    %find distance of each point from this cluster centre
    A = [yi(:),xi(:)];
    B = A; B(:,1) = C(iCluster,1); B(:,2) = C(iCluster,2);
    dx = reshape(nph_haversine(A,B),size(xi));
    clear A B
    
    %discard any points too far from a cluster centre
    dx(dx > Settings.MaxDist) = 1e6;
    
 
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
  disp('Clusters mapped')
  
  clear iCluster Delta dx  C Closer  CD xi yi
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%
  %We now know which point on the map goes with which cluster
  %now find a median for each cluster and map it
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  textprogressbar('Gridding data ')
  for iVar=1:1:numel(Settings.Vars);
    textprogressbar(iVar./numel(Settings.Vars).*100)
    
    %get data
    V = Data.(Settings.Vars{iVar});
    
    %loop over clusters and find median for each
    ClusterMedians = NaN(max(CC(:)),1);
    for iCluster = 1:1:max(CC(:))
      ClusterMedians(iCluster) = nanmedian(V(Data.Cluster == iCluster));
    end; clear iCluster V
    
    %put the data onto the map
    R = Results.(Settings.Vars{iVar});
    Ri = R(iMonth,:,:);    
    for iCluster = 1:1:max(CC(:))
      ThisCluster = find(CC == iCluster);
      Ri(ThisCluster) = ClusterMedians(iCluster);
    end
    R(iMonth,:,:) = Ri;
    Results.(Settings.Vars{iVar}) = R;
    
    
    clear R Ri iCluster ThisCluster
  end; clear iVar  
  
  %also store cluster map
  R = Results.Cid;
  R(iMonth,:,:) = CC;
  Results.Cid = R;
  clear xi yi zz R 
  
  %and number of points used
  R = Results.N;
  Ri = R(iMonth,:,:); Ri(:) = 0;
  for iCluster = 1:1:max(CC(:))
    ThisCluster = find(CC == iCluster);
    N = numel(find(Data.Cluster == iCluster));
    Ri(ThisCluster) = Ri(ThisCluster)+N;
  end
  R(iMonth,:,:) = Ri;
  Results.N = R;
  clear R Ri iCluster  CC ThisCluster
  
  textprogressbar('!')
  
  
  
  %done. save
save(Settings.OutFile,'Results','Settings')
end; clear iMonth

%done. save
save(Settings.OutFile,'Results','Settings')
