clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate monthly maps of IAGOS postprocessed data, composited over all years
%
%forked from parent routine to statistically bootstrap within each cluster
%and work only over a restricted geographic range if desired
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/05/14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%file handling
Settings.DataDir = [LocalDataDir,'/corwin/IAGOS_st/'];
Settings.OutFile = 'v3_test.mat';

%gridding - this is the final map gridding, based on where the clusters
%fall, and does not represent the area averaged over.
Settings.Lon = -180:.5:180;
Settings.Lat = -90:.5:90;

%select pressure ranges
Settings.PrsRanges = [-100,-25;-25,25;25,100]; %hPa

%variables
Settings.Vars = {'STT_A','STT_k','U','V','T','TropPres','Prs'}; %'TropPres','Prs' are required for some internal logic to work, but can be anywhere in the order

%number of clusters
Settings.NClusters = 2000;

%maximum distance of a point from cluster centre (km)
Settings.MaxDist = 500;

%minimum number of points in a cluster to be analysed
Settings.MinPoints = 100;

%bootstrapping parameters
%first is number of samples per strap, second is number of straps
Settings.NSamples = 100;
Settings.Straps   = 100;

%years to use
Settings.Years = 1994:1:2020;

%statistics to compute for each cluster
%numberical values are percentiles, text entries are specific stats
Settings.Stats = {'median','gini','mean','stdev'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create results arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = NaN(12,numel(Settings.Lon),numel(Settings.Lat),numel(Settings.Stats),size(Settings.PrsRanges,1));
Results = struct();
for iVar=1:1:numel(Settings.Vars)
  Results.(Settings.Vars{iVar}) = X;
end;
clear iVar X

Results.Cid = NaN(12,numel(Settings.Lon),numel(Settings.Lat),size(Settings.PrsRanges,1)); %cluster id
Results.N   = Results.Cid ; %number of points actually used in cluster

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
    if mod(iDay,20) == 1; textprogressbar(iDay./numel(Days).*100); end
    
    %find file for this day
    ThisDayFile = wildcardsearch(Settings.DataDir,['*',num2str(Days(iDay)),'*v3*']);
    if numel(ThisDayFile) == 0; clear ThisDayFile; continue; end
    
    %load data
    ThisDayData = load(ThisDayFile{1});
    clear ThisDayFile
    
   
    %discard data outside our geographic region
    InLatRange = inrange(ThisDayData.Results.Lat,[min(Settings.Lat),max(Settings.Lat)]);
    InLonRange = inrange(ThisDayData.Results.Lon,[min(Settings.Lon),max(Settings.Lon)]);
    
    %merge the three criteria
    Good = intersect(InLatRange,      InLonRange);
    
    %store what we need
    for iVar=1:1:numel(Vars);
      V = ThisDayData.Results.(Vars{iVar});
      Data.(Vars{iVar}) = cat(1,Data.(Vars{iVar}),flatten(V(Good)));
    end


    %and done!
    clear ThisDayData iVar V InPrsRange Good InLonRange InLatRange
    
  end
  textprogressbar(100);textprogressbar('!')
  
  
  %my above algorithm leaves some zeros at the end where a longer time
  %series was added to a smaller one. these will all have [lat,lon] = [0,0]. 
  %remove these. We'll also lose any real data at *exactly* [0,0], but this
  %is pretty rare due to the precision of the data (it might never occur.)
  
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
  %% split by tropopause-relative pressure range and then loop over them
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %retain the compelte input data in a special variable
  MASTERDATA = Data; clear Data;
  
  for iLayer = 1:1:size(Settings.PrsRanges,1)
    
    disp(['======== LAYER ',num2str(iLayer),' of ',num2str(size(Settings.PrsRanges,1)),' ============'])
    
    %copy across all the data
    Data = MASTERDATA;
    
    %find the data in this layer
    dTP = Data.Prs - Data.TropPres;
    InRange = find(dTP >= min(Settings.PrsRanges(iLayer,:))...
                 & dTP <= max(Settings.PrsRanges(iLayer,:)));
    Vars = fieldnames(Data);
    for iVar=1:1:numel(Vars)
      Var = Data.(Vars{iVar});
      Data.(Vars{iVar}) = Var(InRange);
    end
    clear Vars iVar Var dTp InRange

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% OK. we have the data lists. Now we need to make the maps
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %k-means cluster the data spatially
    %do this in lat and lon coordinates, as spatial are just as bad due to
    %the large area covered
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    disp('Computing clusters')
    
    %we have far too many points to be tractable in computer time
    %but we know that they're very heavily clumped, since aircraft follow
    %fixed routes. So, let's treat anything in a .2x.2 degree box as
    %'saturated' and as a single point for the cluster generation process
    Lat = round(Data.Lat./.2);
    Lon = round(Data.Lon./.2);
    Combo = Lat + 1000.*Lon;
    [~,idx] = unique(Combo);
    clear Combo Lat Lon
    

    [~,C] = kmeans([Data.Lat(idx),Data.Lon(idx)], ...
      Settings.NClusters, ...
      'Distance','sqeuclidean',...
      'MaxIter',1000, ...
      'Replicates',10, ...
      'Display','off');
     
    %now, find where all the un-clustered points go
    [dx,idx] = pdist2(C,[Data.Lat,Data.Lon],'euclidean','Smallest',1);
    
    %remove points too far away to be fairly included in a cluster
    idx(dx > Settings.MaxDist) = [];
    %store the cluster ids against the data
    Data.Cluster = Data.Lon .* NaN;
    for iCluster = 1:1:size(C,1);
      Data.Cluster(idx == iCluster) = iCluster;
    end; clear iCluster
    
    
    clear idx dx
    
    
    disp('Clusters assigned')
    
    %% assign all geographic points to their nearest cluster centre
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %create results arrays
    CC = NaN( numel(Settings.Lon),numel(Settings.Lat));
    CD = ones(numel(Settings.Lon),numel(Settings.Lat)).*1e7; %1e7km is a fill value that should be larger than any real distance
    
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
      dx(dx > Settings.MaxDist) = 1e7;
      
      
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
    
    textprogressbar('Bootstrapping and gridding data ')
    for iVar=1:1:numel(Settings.Vars);
      textprogressbar((iVar-1)./(2+numel(Settings.Vars)).*100)
      
      %get data
      V = Data.(Settings.Vars{iVar});
      
      %loop over clusters and find requested stats for each
      ClusterStats = NaN(max(CC(:)),numel(Settings.Stats));
      for iCluster = 1:1:max(CC(:))
        
        %pull out the data we have in this cluster, and check it's enough points
        ThisCluster = V(Data.Cluster == iCluster);
        
        %also remove nans. This means we can use non-nansafe routines, which
        %will speed things up a lot at this volume of data
        ThisCluster = ThisCluster(~isnan(ThisCluster));
        
        if numel(ThisCluster) < Settings.MinPoints; continue; end
        
        %bootstrapping time. sample the data first.
        r = randi(numel(ThisCluster),Settings.NSamples,Settings.Straps);
        ThisCluster = ThisCluster(r);
        
        for iStat= 1:1:numel(Settings.Stats)
          
          %declare array for the bootstraps
          Straps = NaN(Settings.Straps,1);
          
          %compute the properties
          switch Settings.Stats{iStat}
            case 'gini';
              %only defined for positive data. It's also really slow. So only use for GW amplitudes.
              if strfind(Settings.Vars{iVar},'ST') & strfind(Settings.Vars{iVar},'_A')
                Straps= ginicoeff(ThisCluster,1);
              end
            case 'mean';
              Straps= mean(ThisCluster,1);
            case 'median';
              Straps= median(ThisCluster,1);
            case 'stdev';
              Straps = std(ThisCluster,[],1);
            otherwise
              if isnumeric(Settings.Stats{iStat});
                Straps = prctile(ThisCluster, Settings.Stats{iStat},1);
              else
                disp('Statistic not specified')
              end
          end
          
          %take the median
          ClusterStats(iCluster,iStat) = nanmedian(Straps);
        end
      end; clear iCluster V Stat ThisCluster iStat idx iStrap
      
      %put the data onto the map
      R = Results.(Settings.Vars{iVar});
      for iStat=1:1:numel(Settings.Stats)
        Ri = R(iMonth,:,:,iStat);
        for iCluster = 1:1:max(CC(:))
          Ri(CC == iCluster) = ClusterStats(iCluster,iStat);
        end
        R(iMonth,:,:,iStat,iLayer) = Ri;
      end
      Results.(Settings.Vars{iVar}) = R;
      
      
      clear R Ri iCluster ThisCluster ClusterStats
    end; clear iVar
    
    %also store cluster map
    R = Results.Cid;
    R(iMonth,:,:,iLayer) = CC;
    Results.Cid = R;
    clear xi yi zz R
    
    %and number of points used
    R = Results.N;
    Ri = R(iMonth,:,:,iLayer); Ri(:) = 0;
    for iCluster = 1:1:max(CC(:))
      ThisCluster = find(CC == iCluster);
      N = numel(find(Data.Cluster == iCluster));
      Ri(ThisCluster) = Ri(ThisCluster)+N;
    end
    R(iMonth,:,:,iLayer) = Ri;
    Results.N = R;
    clear R Ri iCluster  CC ThisCluster
    
    textprogressbar(100);textprogressbar('!')
    
  
  end; clear iLayer
  clear MASTERDATA
  
  %done. save
  save(Settings.OutFile,'Results','Settings','-v7.3')
end; clear iMonth

%done. save
save(Settings.OutFile,'Results','Settings','-v7.3')
