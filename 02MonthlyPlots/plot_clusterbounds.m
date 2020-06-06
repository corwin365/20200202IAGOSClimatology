clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%plot maps of IAGOS properties (including GWs) from prepared data
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/JUN/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%file handling
Settings.Mode = 'h';
Settings.Layers = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load and merge quarterly files
Store = struct();
for iQ=1:1:4;
  
  
  switch iQ; 
    case 1; Q = 'djf'; 
    case 2; Q = 'mam'; 
    case 3; Q = 'jja'; 
    case 4; Q = 'son'; 
  end
  
  Layers = struct();
  for iLayer = 1:1:numel(Settings.Layers)
    File = ['out/',Settings.Mode,'_',Q,'_','b',num2str(Settings.Layers(iLayer)),'.mat'];
    File = load(File);
    if iLayer == 1;
      Layers = File.Results;
    else
      Layers = cat_struct(Layers,File.Results,3);
    end
  end

  if iQ == 1; 
    Store = Layers;
    Meta = File.Settings;
  else
    Store = cat_struct(Store,Layers,5);
  end
  
  
end; clear iQ Q File iLayey Layers

%pull out selected data
Dat2 = struct();
Dat2.Cid = Store.Cid;
Data = Dat2; clear Dat2 iVar Store 

if strcmp(Settings.Mode,'g')
  %need to shift grid-mode data slightly to correctly align geographically
  Meta.Grid.Lat = Meta.Grid.Lat + mean(diff(Meta.ClusterParams.G.Lat))./2;
  Meta.Grid.Lon = Meta.Grid.Lon + mean(diff(Meta.ClusterParams.G.Lon))./2;  
end


disp('Data loaded')



%prepare figure
clf
set(gcf,'color','w')
subplot = @(m,n,p) subtightplot (m, n, p, 0.08, [0.05 0.05], [0.05 0.05]);


% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %% plot 1 - sample map of clusters
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % % subplot(2,2,1)
% % % 
% % % %get data
% % % iQuarter = 1;
% % % ClusterMap = squeeze(Data.Cid(:,:,iQuarter))';
% % % ClusterMap(isnan(ClusterMap)) = 0;
% % % 
% % % %produce map
% % % m_proj('robinson','lat',[-45,85],'lon',[-180 180])
% % % m_coast('patch',[178,255,102]./255,'edgecolor','none'); %coasts at the bottom
% % % hold on
% % % 
% % % %identify all discrete regions in the dataset. these are the individual clusters.
% % % %then plot each as a patch object. also retain areas, for histogram
% % % [xi,yi] = meshgrid(Meta.Grid.Lon,Meta.Grid.Lat);
% % % Area = repmat(112.*mean(diff(Meta.Grid.Lon)) .* 112.*cosd(Meta.Grid.Lat).*mean(diff(Meta.Grid.Lon)),numel(Meta.Grid.Lon),1)';
% % % 
% % % %plot the no-coverage zone
% % % Mask = zeros(size(ClusterMap));
% % % Mask(ClusterMap == 0) = 1; Mask(Mask == 0) = NaN;
% % % colormap([1,1,1;0,0,0])
% % % m_contourf(xi,yi,Mask,1,'edgecolor','none')
% % % 
% % % %plot cluster edges
% % % AreaStore = NaN(max(Data.Cid(:)),1);
% % % textprogressbar('Drawing clusters ')
% % % warning off
% % % for iCluster = 1:1:max(Data.Cid(:));
% % %   if mod(iCluster, 100);textprogressbar(iCluster./max(Data.Cid(:)).*100); end
% % % 
% % %   %identify cluster
% % %   Mask = zeros(size(ClusterMap));
% % %   ThisCluster = find(ClusterMap == iCluster);
% % %   Mask(ThisCluster) = 1;
% % % 
% % %   %find area and save
% % %   AreaStore(iCluster) = nansum(Mask.*Area,[1,2]);
% % %   
% % %   
% % % 
% % %   %plot cluster
% % %   m_contour(xi,yi,Mask,1,'edgecolor','k','linewi',0.25);
% % % 
% % %  
% % % end
% % % textprogressbar('!')
% % % disp('m_grid-ding, this may take a considerable time')
% % % m_grid('fontsize',14)
% % % warning on
% % % 
% % % drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot 2 - sample map of clusters, good region only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,3)

%get data
iQuarter = 1;
ClusterMap = squeeze(Data.Cid(:,:,iQuarter))';
ClusterMap(isnan(ClusterMap)) = 0;

%produce map
m_proj('lambert','lat',[20,80],'lon',[-130 140])
m_coast('patch',[178,255,102]./255,'edgecolor','none'); %coasts at the bottom
hold on

%identify all discrete regions in the dataset. these are the individual clusters.
%then plot each as a patch object. also retain areas, for histogram
[xi,yi] = meshgrid(Meta.Grid.Lon,Meta.Grid.Lat);
Area = repmat(112.*mean(diff(Meta.Grid.Lon)) .* 112.*cosd(Meta.Grid.Lat).*mean(diff(Meta.Grid.Lon)),numel(Meta.Grid.Lon),1)';

%plot the no-coverage zone
Mask = zeros(size(ClusterMap));
Mask(ClusterMap == 0) = 1; Mask(Mask == 0) = NaN;
colormap([1,1,1;0,0,0])
m_contourf(xi,yi,Mask,1,'edgecolor','none')

%plot cluster edges
AreaStore = NaN(max(Data.Cid(:)),1);
textprogressbar('Drawing clusters ')
warning off
for iCluster = 1:1:max(Data.Cid(:));
  if mod(iCluster, 100);textprogressbar(iCluster./max(Data.Cid(:)).*100); end

  %identify cluster
  Mask = zeros(size(ClusterMap));
  ThisCluster = find(ClusterMap == iCluster);
  Mask(ThisCluster) = 1;

  %find area and save
  AreaStore(iCluster) = nansum(Mask.*Area,[1,2]);
  
  

  %plot cluster
  m_contour(xi,yi,Mask,1,'edgecolor','k','linewi',0.25);

 
end
textprogressbar('!')
disp('m_grid-ding, this may take a considerable time')
m_grid('fontsize',14)
warning on
drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot 3 - histogram of cluster area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,2)
cla

%generate histogram
x = linspace(0,3.5,100);
histo = hist(AreaStore./1e5,x);

%plot bars
hold on
for iX=1:1:numel(x)-1
  pos = [x(iX)-mean(diff(x)), 0, x(iX+1)-x(iX), histo(iX)];
  rectangle('position',pos,'facecolor',[1,1,1].*0.7)
%   drawnow
end

% plot(x2,y2)
axis([0 3.5 0 max(histo)*1.05])
xlabel('Cluster Area [x10^5 km^2]'); ylabel('Number of clusters')
hold on

%add a line showing the area of fixed-size boxes
for iGrid = 1:1:6;
  x = 112.*112.*sind(45).*iGrid.*iGrid./1e5;
  plot([1,1].*x,[0,max(histo)*1.1],'k:','linewi',2,'clipping','off')
  text(x+0.02,max(histo).*1.1,[num2str(iGrid),'^\circ'],'clipping','off')
end


box on
drawnow
