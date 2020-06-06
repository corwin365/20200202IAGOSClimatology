clearvars -except SMOOTHSIZE

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

%smoothing of final plot
Settings.SmoothSize = [1,1].*1;%.*SMOOTHSIZE;

%plot rows. one row for each combination of the below
  %%strings, as a cell struct
Settings.Vars    = {'Cid'};%'U','T'};
  %%indices in the order specified in file Settings struct 
Settings.Layers = 1;%7%[7,5]; 
Settings.Stats  = 3;%,4];

Settings.Log = 0;


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
for iVar=1:1:numel(Settings.Vars);
  Dat2.(Settings.Vars{iVar}) = Store.(Settings.Vars{iVar});
end
Data = Dat2; clear Dat2 iVar Store 

if strcmp(Settings.Mode,'g')
  %need to shift grid-mode data slightly to correctly align geographically
  Meta.Grid.Lat = Meta.Grid.Lat + mean(diff(Meta.ClusterParams.G.Lat))./2;
  Meta.Grid.Lon = Meta.Grid.Lon + mean(diff(Meta.ClusterParams.G.Lon))./2;  
end


disp('Data loaded')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%prepare figure
clf
set(gcf,'color','w')
subplot = @(m,n,p) subtightplot (m, n, p, [0.05, 0.02], [0.05 0.05], [0.05 0.05]);


%generate plotting combinations we want to produce
Combos = {};
k = 0;
for iVar=1:1:numel(Settings.Vars)
  for iLayer = 1:1:numel(Settings.Layers)
    for iStat = 1:1:numel(Settings.Stats)
      k = k+1;
      if strcmp('N',  Settings.Vars{iVar}) ...
       | strcmp('Cid',Settings.Vars{iVar}); Combos{k} = {Settings.Vars{iVar},iLayer,1};
      else;Combos{k} = {Settings.Vars{iVar},iLayer,iStat};end
    end
  end
end
clear iVar iLayer iStat k

%begin plotting
k = 0;
for iCombo = 1:1:numel(Combos)
  
  %get data
  Combo = Combos{iCombo};
  ComboData = Data.(Combo{1});
  ComboData = squeeze(ComboData(:,:,Combo{2},Combo{3},:));

  if Settings.Log == 1; ComboData = log10(ComboData); end
  
  
  ColourRange = prctile(ComboData(:),[2.5,97.5]);
  
  for iQuarter = 1:1:4;
    
    %generate panel and map
    k = k+1;
    subplot(2,2,k);
%     subplot(numel(Combos),4,k)
    m_proj('lambert','lat',[20,85],'lon',[-130 140])
%     m_proj('lambert','lat',[30,80],'lon',[-80 40])
    
    %overinterpolate the data so that pcolor/contourf don't hide very thin tendrils
    %using nearest-neighbour interpolant, so no new data introduced
    if ~strcmp(Combo{1},'Cid'); 
      x = Meta.Grid.Lon;
      y = Meta.Grid.Lat;
      [xi,yi] = meshgrid(linspace(min(x),max(x),numel(x).*3), ...
                         linspace(min(y),max(y),numel(y).*3));
      ToPlot = interp2(x,y,squeeze(ComboData(:,:,iQuarter))',xi,yi);
      
      %smooth
      Bad = find(isnan(ToPlot)); ToPlot = inpaint_nans(ToPlot);
      ToPlot = smoothn(ToPlot,Settings.SmoothSize.*5,'gauss',Settings.SmoothSize.*3);
      ToPlot(Bad) = NaN;
    else
      xi = Meta.Grid.Lon; yi = Meta.Grid.Lat; ToPlot = squeeze(ComboData(:,:,iQuarter))'; 
    end
    

    
    
    %plot data
    m_pcolor(xi,yi,ToPlot)
    
    
    
    %tidy up map
    caxis(ColourRange)
    colormap(cbrew('Accent',10))
    
%     colorbar('southoutside')
%     colormap(cbrew('Greens',16));
%     redyellowblue16
    m_coast('color','k');
    m_grid('fontsize',10);
    
    
    %done!
    drawnow
  end
  
  
  
  
end
