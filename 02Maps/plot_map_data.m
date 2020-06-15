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

%smoothing of final plot
Settings.SmoothSize = [1,1].*5; %FWHM of Gaussian smoother

% % %plot rows. one row for each combination of the below
% %   %%strings, as a cell struct
% % Settings.Vars    = {'STT_A','STT_k';};
% %   %%indices in the order specified in file Settings struct 
% % Settings.Layers = 23;
% % Settings.Stats  = [3,4];

%plot rows. one row for each combination of the below
  %%strings, as a cell struct
Settings.Vars    = {'STT_A';};
  %%indices in the order specified in file Settings struct 
Settings.Layers = 22:23;
Settings.Stats  = [1];

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
subplot = @(m,n,p) subtightplot (m, n, p, [0.01, 0.019], [0.05 0.1], [0.15 0.1]);
Letters = 'abcdefghijklmnopqrstuvwxyz';

%generate plotting combinations we want to produce
Combos = {};
k = 0;
for iVar=1:1:numel(Settings.Vars)
  for iLayer = 1:1:numel(Settings.Layers)
    for iStat = 1:1:numel(Settings.Stats)
      k = k+1;
      if k > 3; continue; end
      if strcmp('N',  Settings.Vars{iVar}) ...
       | strcmp('Cid',Settings.Vars{iVar}); Combos{k} = {Settings.Vars{iVar},iLayer,1};
      else;Combos{k} = {Settings.Vars{iVar},iLayer,Settings.Stats(iStat)};end
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
%     subplot(2,2,k);
    hp = subplot(numel(Combos),4,k);
    m_proj('lambert','lat',[25,80],'lon',[-130 150])
%     m_proj('lambert','lat',[30,80],'lon',[-80 40])
    
    %black fill for areas with no data
    for Lon = -180:1:170
      m_patch([0 10 10 0 0]+Lon, ...
              [-90 -90 90 90 -90],....
              [1,1,1].*0.6,'edgecolor','none')
      hold on
    end



    %get data
    xi = Meta.Grid.Lon; yi = Meta.Grid.Lat; ToPlot = squeeze(ComboData(:,:,iQuarter))';
    
    %smooth data
    Bad = find(isnan(ToPlot)); ToPlot = inpaint_nans(ToPlot);
    ToPlot = smoothn(ToPlot,Settings.SmoothSize.*3,'gauss',Settings.SmoothSize./2.355);
    ToPlot(Bad) = NaN;


    %choose colours
    switch Combo{3}
      case {1,3}
        switch Combo{1}
          case 'T';     colormap(hp,cbrew('RdBu',16));           ColourRange = [210,240];  Units = 'Temperature [K]';
          case 'U';     colormap(hp,cbrew('nph_BlueOrange',16)); ColourRange = [-50,50];   Units = 'Zonal Wind [K]';
          case 'STT_A'; colormap(hp,cbrew('RdYlBu',16));         ColourRange = [0.35 1.2]; Units = 'Amplitude [K]';
          case 'STT_k'; colormap(hp,cbrew('Blues',16));           ColourRange = [0 1000]; Units = 'Wavelength [km]'; ToPlot = 1./ToPlot;
        end
      case 4; %only computed upstream for STT_A
         colormap(hp,cbrew('Reds',16)); 
         ColourRange = [0.15 0.45]; 
         Units = 'Gini Coefficient';
    end
    
    
    %change colour scale for gini coefficient
    
    
    %plot data
    m_pcolor(xi,yi,ToPlot)
    
        
    

    
    %tidy up map
    caxis(ColourRange)
    m_coast('color','k');
    m_grid('fontsize',10,'linestyle','--');

    %label
    m_text(0,90,['(',Letters(k),')'],'fontsize',22,'clipping','off', ...
           'horizontalalignment','center','fontweight','bold')
    

    
         %season labels
  if iCombo  == 1
    switch iQuarter
      case 1; Q = 'djf';
      case 2; Q = 'mam';
      case 3; Q = 'jja';
      case 4; Q = 'son';
    end
    m_text(0,150,upper(Q),'horizontalalignment','center','fontsize',35,'fontweight','bold')
  end
  
  %done!
  drawnow
    
  end
  
  
  %colour table
  x = [0.92 0.02];
  y = [-1,1].*.1 + (1./(numel(Combos)+1))*((numel(Combos)+1)-iCombo) + (0.5-(iCombo-numel(Combos)./2)).*0.03;
  cb = colorbar('position',[x(1),y(1),x(2),y(2)-y(1)]);
  cb.Label.String = Units;
  drawnow


  
end

