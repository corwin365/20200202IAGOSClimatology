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
Settings.DataFile = 'allheights_c3000_s10mil_80_500.mat';;%'allheights_c3000_s10mil.mat';%'allheights_c3000_s10mil.mat';

%variable to plot
Settings.Var = 'T';

%statistic to plot (number in order of input file)
Settings.Stat = 1; %ignored for N and Cid

%smoothing (bins)
Settings.SmoothSize =[1,1].*11;

%colours
Settings.NColours = 16;

%gap filling. maximum number of bins permitted for a fill
%this happens *before* smoothing
%not used as clustering plugs the same bags, but vestigial code kept
Settings.MaxGapSize.Lon = 0;
Settings.MaxGapSize.Lat = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load and prep data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FIGCOUNT = 0;
Letters = 'abcdefghijkl';
%plot settings
clf
set(gcf,'color','w')
subplot = @(m,n,p) subtightplot (m, n, p, 0.04, 0.05, [0.03,0.12]);
  

for LAYER =1;% 1:1:3;

  Settings.Layer = LAYER;
  switch LAYER
    case 1 ; LayerName = 'Lower Stratosphere';
    case 2 ; LayerName = 'Tropopause';
    case 3 ; LayerName = 'Upper Troposphere';
  end
  
  %load file
  Data = load(Settings.DataFile);
  
  %find desired var and stat
  Data.Results = Data.Results.(Settings.Var);
  
  if strcmp(Settings.Var,'N') | strcmp(Settings.Var,'Cid');
    Data.Results = Data.Results(:,:,:,Settings.Layer);    
  else
    Data.Results = Data.Results(:,:,:,Settings.Stat,Settings.Layer);
  end
  
  %special handling
  switch Settings.Var
    case 'Prs';   Data.Results = Data.Results ./ 100;
    case 'STT_k'; Data.Results = 1./Data.Results;
  end
  
  
  
% %   %colour range
% %   switch Settings.Var
% %     case {'U','V'};        CRange = [-1,1].*prctile(abs(Data.Results(:)),97.5);
% %     case {'dW'};           CRange = [-1,1].*prctile(abs(Data.Results(:)),66);
% %     case {'Prs','W'};      CRange = prctile(abs(Data.Results(:)),[2.5,97.5]);
% %     case {'STT_k'};        CRange = [50,300];
% %     case {'T'};            CRange = [210 234];
% %     case {'STT_A'};        CRange = [0,1.5];
% %     otherwise;             CRange = [0,prctile(Data.Results(:),97.5)];
% %   end
  
  CRange = prctile(abs(Data.Results(:)),[2.5,97.5]);
  
% %   %load topography
% %   Topo = topo_etc([-180,179],[-90,89]);
% %   
% %   %downsample topography
% %   lon2 = -180:0.5:180;
% %   lat2 = -90:0.5:90;
% %   [xi,yi] = meshgrid(lon2,lat2);
% %   Topo.elev = interp2(Topo.lons,Topo.lats,Topo.elev,xi,yi);
% %   Topo.lons = xi; Topo.lats = yi; clear xi yi
  
  
  %name that stat!
  if strcmp(Settings.Var,'N') | strcmp(Settings.Var,'Cid')
    StatName = '';
  else
    switch Data.Settings.Stats{Settings.Stat}
      case 'gini';    StatName = 'Gini Coefficient of';
      case 'mean';    StatName = 'Mean';
      case 'stdev';   StatName = 'Standard Deviation of';
      case 'median';  StatName = 'Median';
      otherwise;
        if isnumeric(Data.Settings.Stats{Settings.Stat});
          StatName = [num2str(Data.Settings.Stats{Settings.Stat}),'th %ile of'];
        else;
          StatName = '???';
        end
    end
  end
  
  %name that variable
  switch Settings.Var
    case 'A';    VarName = 'Wave Amplitude';
    case 'k';    VarName = 'Along-flight Wavelength';
    case 'T';    VarName = 'Temperature';
    case 'U';    VarName = 'Zonal Wind Speed';
    case 'V';    VarName = 'Merid. Wind Speed';
    case 'Prs';  VarName = 'Air Pressure';
    otherwise;   VarName = Settings.Var;
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% plot each month
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

  for iQuarter=3%:1:4
    FIGCOUNT = FIGCOUNT+1;
    
    %create subplot
%       subplot(3,4,FIGCOUNT)
%       subplot(2,2,FIGCOUNT)
%     subplot(1,4,FIGCOUNT)
    
    %plot settings
    cla
    
    %create map
%     m_proj('stereographic','lat',90,'long',0,'radius',70);
     m_proj('lambert','lon',[-130,170],'lat',[20,80]);
%      m_proj('robinson','lon',[-130,180],'lat',[-40,90]);
    
    %get data
    ToPlot = squeeze(Data.Results(iQuarter,:,:))';
    
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
    
%     %no colours below table
%     ToPlot(ToPlot < min(CRange)) = min(CRange);
    
    %colour levels
    CLevels = linspace(min(CRange),max(CRange),Settings.NColours+1);
    
    
    %plot data
    %   m_contourf(Data.Settings.Lon, ...
    %              Data.Settings.Lat, ...
    %              ToPlot,CLevels,'edgecolor','none');
    m_pcolor(Data.Settings.Lon,Data.Settings.Lat,ToPlot);
    hold on
    shading flat
    
    %colours
    % %   switch Settings.Var;
    % %     case 'A';  colormap(cbrew('YlOrRd',4));
    % %     otherwise; colormap(cbrew('rdYlBu',4));
    % %   end
    colormap(cbrew('RdYlBu',Settings.NColours));
%     colorbar('southoutside')
%     caxis(CRange)
% caxis([0.2 1.5])
    
% %     %plot topography
% %     m_contour(Topo.lons,Topo.lats,Topo.elev,1:2:7,'color',[1,1,1].*0.4,'linewi',0.5);
    
    %tidy
    m_coast('color',[1,1,1].*0.2,'linewi',0.5);
    m_grid('fontsize',14,'ytick',[],'xtick',[-120:60:180]);

    %label
    m_text(-5,90,['(',Letters(FIGCOUNT),')'],'fontsize',30,'clipping','off','fontweight','bold','horizontalalign','center')
    
    drawnow

    
  end
end
  %%
cb = colorbar('eastoutside','Position',[0.90,0.36,0.02,0.3]);
cb.Label.String = [StatName,' ',VarName];
