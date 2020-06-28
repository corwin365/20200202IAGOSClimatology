clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot ERA5 maps  of basic vars which can be compared to ours
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings.Era5File = 'era5_neilmap.mat';
Settings.Era5File = 'era5_250_200.mat';
Settings.MapFile = '../02Maps/out/h_DJF_b1_sgolay900.mat'; %mjust an example to pull the grid from

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Era5 = load(Settings.Era5File);
Map  = load(Settings.MapFile);

%latitudes are backwards in the ERA5 files. Oops, my fault.
Era5.Settings.LatScale = Era5.Settings.LatScale(end:-1:1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sz = size(Map.Results.T);
NewMaps = NaN([3,sz]);
Diff    = NewMaps;
OldMaps = NewMaps;

[~,mm,~] = datevec(Era5.Settings.TimeScale);
for iQ = 1:1:4
  disp(iQ)
  switch iQ
    case 1; InRange = find(mm ==  9 | mm == 10 | mm == 11); Map = load(['../02Maps/out/h_SON_b1_sgolay900.mat']);
    case 2; InRange = find(mm == 12 | mm ==  1 | mm ==  2); Map = load(['../02Maps/out/h_DJF_b1_sgolay900.mat']);
    case 3; InRange = find(mm ==  3 | mm ==  4 | mm ==  5); Map = load(['../02Maps/out/h_MAM_b1_sgolay900.mat']);
    case 4; InRange = find(mm ==  6 | mm ==  7 | mm ==  8); Map = load(['../02Maps/out/h_JJA_b1_sgolay900.mat']);
  end
  
  for iVar=1:1:3;
    
    %get data
    switch iVar
      case 2; Data = Era5.Results.U;  MapO = Map.Results.U(:,:,1,3);
      case 3; Data = Era5.Results.V;  MapO = Map.Results.V(:,:,1,3);
      case 1; Data = Era5.Results.T;  MapO = Map.Results.T(:,:,1,3);
    end
    
    %select time range and take median
    Data = squeeze(nanmedian(Data(InRange,:,:),1));
    
    %interpolate to our output grid
    [xi,yi] = meshgrid(Map.Settings.Grid.Lat,Map.Settings.Grid.Lon);
    Data2 = interp2(Era5.Settings.LatScale,Era5.Settings.LonScale,...
                    Data,xi,yi);
               
    MapO  = smoothn2(MapO,[5,5].*3,'gauss',[5,5]./2.355);
    Data2 = smoothn2(Data2,[5,5].*3,'gauss',[5,5]./2.355);
    
    %store values and difference
    NewMaps(iVar,:,:,1,iQ) = Data2;
    Diff(   iVar,:,:,1,iQ) = MapO - Data2; 
    OldMaps(iVar,:,:,1,iQ) = MapO;
    
  end 
end

Lat = Map.Settings.Grid.Lat;
Lon = Map.Settings.Grid.Lon;
clearvars -except NewMaps Diff Lat Lon OldMaps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iForm = 1:1:3;
  
  figure(iForm); clf
  set(gcf,'color','w')
  subplot = @(m,n,p) subtightplot (m, n, p, [0.01, 0.019], [0.05 0.1], [0.15 0.1]);
  Letters = 'abcdefghijklmnopqrstuvwxyz';
  
  k = 0;
  for iVar=1:1:3
    for iQ=1:1:4;
      
      switch iForm
        case 1; ToPlot = squeeze(OldMaps(iVar,:,:,1,iQ));      
        case 2; ToPlot = squeeze(Diff(   iVar,:,:,1,iQ));    
        case 3; ToPlot = squeeze(NewMaps(iVar,:,:,1,iQ));
      end
          
      k = k+1;
      hp = subplot(3,4,k);
      
      
      m_proj('lambert','lat',[25,80],'lon',[-130 150])
%       m_proj('stereographic','lat',90,'lon',0,'radius',70)
      m_pcolor(Lon,Lat,ToPlot')
      m_coast('color','k');
      m_grid('fontsize',10,'linestyle','--');
      
      switch iVar
        case 1;     colormap(hp,cbrew('RdBu',16));  Units = 'Temperature [K]';           
          if iForm == 1 | iForm == 3; ColourRange = [210,235]; else; ColourRange = [-5,5]; end
        case 2;     colormap(hp,cbrew('nph_BlueOrange',16));  Units = 'Zonal Wind [m/s]';
          if iForm == 1 | iForm == 3; ColourRange = [-50,50]; else; ColourRange = [-10,10]; end 
        case 3;     colormap(hp,cbrew('nph_BlueOrange',16)); Units = 'Merid Wind [m/s]';
          if iForm == 1 | iForm == 3; ColourRange = [-20,20]; else; ColourRange = [-10,10]; end
      end

      caxis(ColourRange);
      
      m_text(0,90,['(',Letters(k),')'],'fontsize',22,'clipping','off', ...
        'horizontalalignment','center','fontweight','bold')
      
      %colour table
      if iQ == 4;
        x = [0.92 0.02];
        y = [-1,1].*.1 + (1./4)*(4-iVar) + (0.5-(iVar-3./2)).*0.03;
        cb = colorbar('position',[x(1),y(1),x(2),y(2)-y(1)]);
        cb.Label.String = Units;
      end
      
      if iVar == 1;
        switch iQ
          case 2; Q = 'djf';
          case 3; Q = 'mam';
          case 4; Q = 'jja';
          case 1; Q = 'son';
        end
        m_text(0,150,upper(Q),'horizontalalignment','center','fontsize',35,'fontweight','bold')
      end
      
      drawnow
      
    end
  end
  
end