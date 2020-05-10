clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot monthly maps of IAGOS postprocessed data, composited over all years
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/05/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%file handling
Settings.DataFile = 'IAGOS_maps_200hPa.mat';

%variable to plot
Settings.Var = 'T';

%smoothing
Settings.SmoothSize = [3,3];

%gap filling. maximum number of bins permitted for a fill
%this happens *before* smoothing
Settings.MaxGapSize.Lon = 3;
Settings.MaxGapSize.Lat = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load and prep data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load file
Data = load(Settings.DataFile);

%find desired var
Data.Results = Data.Results.(Settings.Var);

%special handling
switch Settings.Var
  case 'Prs'; Data.Results = Data.Results ./ 100;
  case 'k';   Data.Results = 1./Data.Results;
end



%colour range
switch Settings.Var
  case {'U','V'};           CRange = [-1,1].*prctile(abs(Data.Results(:)),95);
  case {'T','Prs','k','A'}; CRange = prctile(abs(Data.Results(:)),[5,95]);    
  otherwise;                CRange = [0,prctile(Data.Results(:),95)];
end  

%load topography
Topo = topo_etc([-180,179],[-90,89]);

%downsample topography
lon2 = -180:0.5:180;
lat2 = -90:0.5:90;
[xi,yi] = meshgrid(lon2,lat2);
Topo.elev = interp2(Topo.lons,Topo.lats,Topo.elev,xi,yi);
Topo.lons = xi; Topo.lats = yi; clear xi yi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot each month
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot settings
clf
set(gcf,'color','w')
subplot = @(m,n,p) subtightplot (m, n, p, 0.03, 0.025, 0.025);

for iMonth=1:1:12
  
  %create subplot
  subplot(3,4,iMonth)
  
  %plot settings
  cla
  
  %create map
  m_proj('lambert','lon',[-125,35],'lat',[20,80]);

  %get data
  ToPlot = squeeze(Data.Results(iMonth,:,:))';
  
  %interpolate over small gaps (only have a 1d routine to do this)
  %first, gaps in longitude
  if Settings.MaxGapSize.Lon ~= 0;
    %loop over rows
    for iRow=1:1:numel(Data.Settings.Lat)
      %extract row
      Row = ToPlot(iRow,:);
      %check the row is not empty
      if sum(~isnan(Row)) < 2; continue; end
      %interpolate over small gaps
      Row = interp1gap(Data.Settings.Lon,Row,Data.Settings.Lon,Settings.MaxGapSize.Lon.*mean(diff(Data.Settings.Lon)));
      %and store
      ToPlot(iRow,:) = Row;
    end
  end
  
  %second, gaps in latitude
  if Settings.MaxGapSize.Lat ~= 0;
    %loop over rows
    for iRow=1:1:numel(Data.Settings.Lon)
      %extract row
      Row = ToPlot(:,iRow);
      %check the row is not empty
      if sum(~isnan(Row)) < 2; continue; end
      %interpolate over small gaps
      Row = interp1gap(Data.Settings.Lat,Row,Data.Settings.Lat,Settings.MaxGapSize.Lat.*mean(diff(Data.Settings.Lon)));
      %and store
      ToPlot(:,iRow) = Row;
    end
  end
  
  
  %smooth
  Bad = find(isnan(ToPlot));
  ToPlot = smoothn(inpaint_nans(ToPlot),Settings.SmoothSize);
  ToPlot(Bad) = NaN;
  
  
  %plot data
  m_contourf(Data.Settings.Lon, ...
             Data.Settings.Lat, ...
             ToPlot,linspace(min(CRange),max(CRange),12),'edgecolor','none');
  hold on
  shading flat

  %colours
  switch Settings.Var;
    case 'A';  colormap(cbrew('YlOrRd',12));
    otherwise; colormap(cbrew('rdYlBu',12));
  end
  caxis(CRange)
  
  %plot topography
  m_contour(Topo.lons,Topo.lats,Topo.elev,0:1:7,'color',[1,1,1].*0.3,'linewi',0.5);
  
  %tidy
  m_coast('color',[1,1,1].*0,'linewi',1);
  m_grid('fontsize',8);
  title(datestr(datenum(2000,iMonth,1),'mmmm'))
  drawnow
  
  
end
%%
colorbar('southoutside','Position',[0.42,0.645,0.16,0.03])
sgtitle(Settings.Var,'fontsize',24)