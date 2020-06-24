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
Settings.SmoothSize = [1,1].*3; %FWHM of Gaussian smoother

%plot rows. one row for each combination of the below
  %%strings, as a cell struct
Settings.Vars    = {'STT_A','U','STT_k','TropPres','V'};
  %%indices in the order specified in file Settings struct 
Settings.Layers = 1;
Settings.Stats  = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load and merge quarterly files
Store = struct();
for iQ=1:1:4;
  
  
  switch iQ; 
    case 2; Q = 'djf'; 
    case 3; Q = 'mam'; 
    case 4; Q = 'jja'; 
    case 1; Q = 'son'; 
  end
  
  Layers = struct();
  for iLayer = 1:1:numel(Settings.Layers)
    File = ['out/',Settings.Mode,'_',Q,'_','b',num2str(Settings.Layers(iLayer)),'_sgolay900.mat'];
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

if isfield(Data,'STT_k'); Data.STT_k =1./Data.STT_k; end


disp('Data loaded')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%prepare figure
clf
set(gcf,'color','w')
subplot = @(m,n,p) subtightplot (m, n, p, [0.01, 0.019], [0.05 0.1], [0.15 0.1]);
Letters = 'abcdefghijklmnopqrstuvwxyz';

%begin plotting
k = 0;


for iQuarter = 1:1:4;
  
  %generate panel and map
  hp = subplot(2,2,iQuarter);
  m_proj('lambert','lat',[25,80],'lon',[-130 150])
  
  %black fill for areas with no data
  for Lon = -180:1:170
    m_patch([0 10 10 0 0]+Lon, ...
      [-90 -90 90 90 -90],....
      [1,1,1].*0.6,'edgecolor','none')
    hold on
  end
  
  %get data
  xi = Meta.Grid.Lon; yi = Meta.Grid.Lat;
  ToPlot = squeeze(Data.STT_A(:,:,1,Settings.Stats,iQuarter))';
  Wind   = squeeze(Data.U(    :,:,1,Settings.Stats,iQuarter))'; 
  Lh     = squeeze(Data.STT_k(:,:,1,Settings.Stats,iQuarter))';  
  
  %smooth data
  ToPlot = smoothn2(ToPlot,Settings.SmoothSize.*3,'gauss',Settings.SmoothSize./2.355);
  Wind   = smoothn2(Wind,  Settings.SmoothSize.*3,'gauss',Settings.SmoothSize./2.355);
  Lh     = smoothn2(Lh,    Settings.SmoothSize.*3,'gauss',Settings.SmoothSize./2.355);
  
  %choose colours
  colormap(hp,cbrew('RdYlBu',16));
%   ColourRange = [0.35 1.2];
% ColourRange = [150 350];
%   Units = 'Amplitude [K]';

  
  %plot data
  ToPlot(Wind < 20 & Lh > 500) = NaN;
  
  m_pcolor(xi,yi,ToPlot); hold on
  colorbar
  
  %overplot wind
%   for U=-50:10:-10
%     [c,h] = m_contour(xi,yi,Wind,[U,U],'k-','linewi',-U./10);
%     clabel(c,h)
%   end  
%   for U=20;%10:10:50
%     [c,h] = m_contour(xi,yi,Wind,[U,U],'k-','linewi',U./10);
%     clabel(c,h)
%   end
  
  %tidy up map
%   caxis(ColourRange)
colorbar
  m_coast('color',[1,1,1].*0.5);
  m_grid('xtick',[],'ytick',[]);%'fontsize',10,'linestyle','--');
  
  %label
  m_text(0,90,['(',Letters(iQuarter),')'],'fontsize',22,'clipping','off', ...
    'horizontalalignment','center','fontweight','bold')
  
  
% %   
% %   %season labels
% %   switch iQuarter
% %     case 2; Q = 'djf';
% %     case 3; Q = 'mam';
% %     case 4; Q = 'jja';
% %     case 1; Q = 'son';
% %   end
% %   m_text(0,150,upper(Q),'horizontalalignment','center','fontsize',35,'fontweight','bold')
  
  %done!
  drawnow
  
end


% % %colour table
% % x = [0.92 0.02];
% % y = [-1,1].*.1 + (1./(numel(Combos)+1))*((numel(Combos)+1)-iCombo) + (0.5-(iCombo-numel(Combos)./2)).*0.03;
% % cb = colorbar('position',[x(1),y(1),x(2),y(2)-y(1)]);
% % cb.Label.String = Units;
% % drawnow
