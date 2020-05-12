clearvars


disp('Does not yet apply alpha')
stop
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
Settings.DataFile = 'clustertest.mat';%'IAGOS_maps_200hPa_median.mat';

%variable to plot
Settings.Var = 'A';

%smoothing (bins)
Settings.SmoothSize =[1,1].*1;

%colours
Settings.NColours = 16;

%gap filling. maximum number of bins permitted for a fill
%this happens *before* smoothing
Settings.MaxGapSize.Lon = 5;
Settings.MaxGapSize.Lat = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load and prep data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load file
Data = load(Settings.DataFile);

%compute the measurements per km^2 for each cluster
%then assign that to each gridbox
[xi,yi] = meshgrid(Data.Settings.Lon,Data.Settings.Lat);
A = 112 .* 112. *cosd(yi);  %km^2 per gridbox


Mpkm = NaN.*Data.Results.N;
for iMonth=1:1:12;
  
  %pull out clister ids for month
  Cid = squeeze(Data.Results.Cid(iMonth,:,:));
  N   = squeeze(Data.Results.N(  iMonth,:,:));
  if nansum(Cid == 0); continue; end
  
  ThisMonth = NaN.*Cid;
  for Cluster = 1:1:nanmax(Cid(:))
    
   
    %find all gridboxes with this cluster ID
    ThisCid = find(Cid == Cluster);
    if numel(ThisCid) == 0; continue; end
    
    %sum their areas
    Area = sum(A(ThisCid));
    
    %and the number of contibuting measurements
    Ntm = mean(N(ThisCid)); %they're all the same number, so mean is fine
    
    %hence, measurements per unit area
    ThisMonth(ThisCid) = Ntm./Area;
  end
  %and store
  Mpkm(iMonth,:,:) = ThisMonth;
end
clear xi yi A iMonth Cid N ThisMonth ThisCid Cluster Area Ntm



%find desired var
Data.Results = Data.Results.(Settings.Var);

%special handling
switch Settings.Var
  case 'Prs'; Data.Results = Data.Results ./ 100;
  case 'k';   Data.Results = 1./Data.Results;
end



%colour range
switch Settings.Var
  case {'U','V'};           CRange = [-1,1].*prctile(abs(Data.Results(:)),97.5);
  case {'T','Prs','k','A'}; CRange = prctile(abs(Data.Results(:)),[2.5,97.5]);    
  otherwise;                CRange = [0,prctile(Data.Results(:),97.5)];
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
  m_proj('lambert','lon',[-130,130],'lat',[20,77]);
%   m_proj('stereographic','lat',90,'long',0,'radius',90);

  %get data
  ToPlot = squeeze(Data.Results(iMonth,:,:))';
  
  %get alpha, generated from the measurements per km^2
  Alpha = squeeze(Mpkm(iMonth,:,:))';
  
  %more than 0.05 measurements per km^2 is opaque, then scales down to zero
  %below that value
  Alpha = Alpha ./ 0.05;
  Alpha(Alpha > 1) = 1;
  
  
  CRange = prctile(abs(Mpkm(:)),[2.5,97.5]);    
  
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
  
  %no colours below table
  ToPlot(ToPlot < min(CRange)) = min(CRange);
  
 
  
  %plot data
% %   m_contourf(Data.Settings.Lon, ...
% %              Data.Settings.Lat, ...
% %              ToPlot,linspace(min(CRange),max(CRange),Settings.NColours),'edgecolor','none');
  m_pcolor(Data.Settings.Lon,Data.Settings.Lat,ToPlot);
  hold on
  shading flat

  %colours
% %   switch Settings.Var;
% %     case 'A';  colormap(cbrew('YlOrRd',12));
% %     otherwise; colormap(cbrew('rdYlBu',12));
% %   end
  colormap(cbrew('RdYlBu',Settings.NColours));
%   caxis(CRange)
  
  %plot topography
  m_contour(Topo.lons,Topo.lats,Topo.elev,1:2:7,'color',[1,1,1].*0.4,'linewi',0.5);
  
  %tidy
  m_coast('color',[1,1,1].*0.2,'linewi',0.5);
  m_grid('fontsize',8);
  title(datestr(datenum(2000,iMonth,1),'mmmm'))
  drawnow
  
  
end
%
colorbar('southoutside','Position',[0.42,0.66,0.16,0.03])
sgtitle(Settings.Var,'fontsize',24)