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

%height regions (relative to tropopause, defined as acp-17-12495-2017)
Settings.Region.LMS = [-100,-25];  %lowermost stratosphere
Settings.Region.TPL = [-25,25];   %tropopause layer
Settings.Region.UTr = [25,100];   %upper troposphere
Settings.Region.All = [-999,999]; %all data

%number of colours
Settings.NColours = 32;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load data and split by region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load data
Data = load(Settings.InFile);

%define regions
idx.LMS = find(Data.Settings.Grid.dTP >= min(Settings.Region.LMS) & Data.Settings.Grid.dTP <= max(Settings.Region.LMS));
idx.TPL = find(Data.Settings.Grid.dTP >= min(Settings.Region.TPL) & Data.Settings.Grid.dTP <= max(Settings.Region.TPL));
idx.UTr = find(Data.Settings.Grid.dTP >= min(Settings.Region.UTr) & Data.Settings.Grid.dTP <= max(Settings.Region.UTr));
idx.All = find(Data.Settings.Grid.dTP >= min(Settings.Region.All) & Data.Settings.Grid.dTP <= max(Settings.Region.All));

%hence, make three calendars
Map.LMS = nansum(Data.Results.Map.dTP(:,:,idx.LMS),3);
Map.TPL = nansum(Data.Results.Map.dTP(:,:,idx.TPL),3);
Map.UTr = nansum(Data.Results.Map.dTP(:,:,idx.UTr),3);
Map.All = nansum(Data.Results.Map.dTP(:,:,idx.All),3);

Lon = Data.Settings.Grid.Lon;
Lat = Data.Settings.Grid.Lat;

%divide maps by total
Map.LMS = Map.LMS ./ nansum(Map.LMS(:));
Map.TPL = Map.TPL ./ nansum(Map.TPL(:));
Map.UTr = Map.UTr ./ nansum(Map.UTr(:));
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
    case 1; ThisMap = Map.UTr; title('(e) Upper Troposphere');
    case 2; ThisMap = Map.TPL; title('(f) Tropopause Layer');
    case 3; ThisMap = Map.LMS; title('(g) Lowermost Stratosphere');      
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


