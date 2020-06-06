clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot IAGOS postprocessed metadata - time distribution
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/05/14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%file handling
Settings.InFile = 'metadata_all_v2.mat';

%height regions (absolute pressure)
Settings.Region.Bot = [1000,225];
Settings.Region.Mid = [225,205];
Settings.Region.Top = [205,0];
Settings.Region.All = [1000,0];

%number of colours
Settings.NColours = 32;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load data and split by region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load data
Data = load(Settings.InFile);


%define regions
idx.Bot = find(Data.Settings.Grid.Prs >= min(Settings.Region.Bot) & Data.Settings.Grid.Prs <= max(Settings.Region.Bot));
idx.Mid = find(Data.Settings.Grid.Prs >= min(Settings.Region.Mid) & Data.Settings.Grid.Prs <= max(Settings.Region.Mid));
idx.Top = find(Data.Settings.Grid.Prs >= min(Settings.Region.Top) & Data.Settings.Grid.Prs <= max(Settings.Region.Top));
idx.All = find(Data.Settings.Grid.Prs >= min(Settings.Region.All) & Data.Settings.Grid.Prs <= max(Settings.Region.All));

%hence, make three calendars
Map.Bot = nansum(Data.Results.Map.Prs(:,:,idx.Bot),3);
Map.Mid = nansum(Data.Results.Map.Prs(:,:,idx.Mid),3);
Map.Top = nansum(Data.Results.Map.Prs(:,:,idx.Top),3);
Map.All = nansum(Data.Results.Map.Prs(:,:,idx.All),3);


Lon = Data.Settings.Grid.Lon;
Lat = Data.Settings.Grid.Lat;

%divide maps by total
Map.Bot = Map.Bot ./ nansum(Map.Bot(:));
Map.Mid = Map.Mid ./ nansum(Map.Mid(:));
Map.Top = Map.Top ./ nansum(Map.Top(:));
Map.All = Map.All ./ nansum(Map.All(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define colours
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Colours = cbrew('RdYlBu',Settings.NColours);
Range = [-3,0.25];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot data. do this by hand with patches, for more control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%prepare figure
clf
set(gcf,'color','w')
subplot = @(m,n,p) subtightplot (m, n, p, 0.028, 0.1, [0.03,0.08]);

%plot data
for iRegion=1:1:4;
  
  %create panel
  subplot(1,4,iRegion)
  m_proj('Robinson','lon',[-180,175],'lat',[-45,85])
  hold on
  
  %get data
  switch iRegion
    case 1; ThisMap = Map.Bot; title('(e) Below 225 hPa');
    case 2; ThisMap = Map.Mid; title('(f) 225 hPa - 205 hPa');
    case 3; ThisMap = Map.Top; title('(g) Above 205 hPa');      
    case 4; ThisMap = Map.All; title('(h) All Data');          
  end
  
  
% %   
  
  %overinterpolate the data with nearest neighbour
  %this is so that we can see all the points, as pcolor shades across
  %the diagonal
  Lon2 = -180:0.1:180;
  Lat2 = -90:0.1:90;
  [Lon2,Lat2] = meshgrid(Lon2,Lat2);
  ThisMap2 = interp2(Lon,Lat,ThisMap',Lon2,Lat2,'nearest')';
  ThisMap2(ThisMap2 == 0) = NaN;
  
  %plot data
  m_pcolor(Lon2,Lat2,log10(ThisMap2'.*100))
  
  %tidy up
  m_coast('color',[1,1,1].*0.5,'linewi',0.5);
  m_grid('fontsize',13,'xtick',[-120,0,120]);
  colormap(Colours)
  caxis(Range)
  
  drawnow
  
 
end


%% done! colourbar
colormap(Colours)
cb = colorbar('position',[0.952,0.4,0.01,0.2]);
caxis(Range)
cb.Label.String = '% of points in gridbox';
set(cb,'ytick',-4:1:0,'yticklabel',{'10^{-4}','10^{-3}','10^{-2}','0.1','1'})


