clearvars


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot seasonal zonal mean-height plots of GW activity overlaid with wind speed
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/06/03
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Settings.InFile = 'zm_longrange.mat';


Settings.Seasons = [12,1,2;3,4,5;6,7,8;9,10,11];

Settings.SmoothSize = [1,1].*3; %in raw data bins
Settings.InterpFactor = [5,5]; 
    
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prepare data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load file
Data = load(Settings.InFile);

%pull out scales and titles
Lat     = Data.Settings.LatScale;
Time    = Data.Settings.TimeScale;
Prs.Abs = Data.Settings.AbsPrsScale;
Prs.Rel = Data.Settings.RelPrsScale;

%generate plot data
Plot.Abs.U     = NaN(4,numel(Lat),numel(Prs.Abs));
Plot.Rel.U     = NaN(4,numel(Lat),numel(Prs.Rel));
Plot.Abs.STT_A = NaN(4,numel(Lat),numel(Prs.Abs));
Plot.Rel.STT_A = NaN(4,numel(Lat),numel(Prs.Rel));
Plot.TP        = NaN(4,numel(Lat));

%months, for use in season identification
[~,mm,~] = datevec(Time);

for iSeason = 1:1:4;
  
  %identify days in season
  ThisSeason = find(ismember(mm,Settings.Seasons(iSeason,:)));
  
  %average wrt time and store
  Plot.Abs.U(  iSeason,:,:) = squeeze(nanmean(Data.Results.Abs(2,ThisSeason,:,:),2));
  Plot.Rel.U(  iSeason,:,:) = squeeze(nanmean(Data.Results.Rel(2,ThisSeason,:,:),2));

  Plot.Abs.STT_A(iSeason,:,:) = squeeze(nanmean(Data.Results.Abs(3,ThisSeason,:,:),2));
  Plot.Rel.STT_A(iSeason,:,:) = squeeze(nanmean(Data.Results.Rel(3,ThisSeason,:,:),2));  
  
  %tropopause height is done separately either way
  Plot.TP(iSeason,:) = nanmean(Data.Results.TP(ThisSeason,:));
  
  
  
end; clear iSeason mm Time ThisSeason 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% special handling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%colour table
Settings.ColourTable = 'nph_RdBuPastel';



%colour range
ColourRange = prctile([Plot.Abs.STT_A(:);Plot.Rel.STT_A(:)],[5,95]);


%rescale smoothing for overinterpolation
Settings.SmoothSize = Settings.InterpFactor.*Settings.SmoothSize;
if ~isodd(Settings.SmoothSize(1)); Settings.SmoothSize(1) = Settings.SmoothSize(1)+1; end
if ~isodd(Settings.SmoothSize(2)); Settings.SmoothSize(2) = Settings.SmoothSize(2)+1; end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf
set(gcf,'color','w')
subplot = @(m,n,p) subtightplot (m, n, p, [0.05, 0.02], [0.1 0.1], [0.05 0.15]);

Letters = 'abcdefgh';
k = 0;

for iBasis=1:1:2;
  
  switch iBasis
    case 1; Data = Plot.Abs.STT_A; Scale = Prs.Abs; YLabel = 'Pressure [hPa]';
    case 2; Data = Plot.Rel.STT_A; Scale = Prs.Rel; YLabel = '\Delta Pressure [hPa]';
  end
      
  switch iBasis
    case 1; AllWind = Plot.Abs.U; Scale = Prs.Abs;
    case 2; AllWind = Plot.Rel.U; Scale = Prs.Rel;
  end

  
  for iQuarter = 1:1:4;
    
    k = k+1; subplot(2,4,k);
    
    %pull data
    ToPlot = squeeze(   Data(iQuarter,:,:))';
    Wind   = squeeze(AllWind(iQuarter,:,:))';
    
    %overinterpolate the data. This is so that pcolor/countourf show
    %data when one vertex is NaN
    x = linspace(  min(Lat),  max(Lat),ceil(range(  Lat)./(mean(diff(  Lat))./Settings.InterpFactor(1))));
    y = linspace(min(Scale),max(Scale),ceil(range(Scale)./(mean(diff(Scale))./Settings.InterpFactor(2))));
    [xi,yi] = meshgrid(x,y);
    ToPlot = interp2(Lat,Scale,ToPlot,xi,yi,'nearest');
    Wind   = interp2(Lat,Scale,Wind,xi,yi,'nearest');
    
    %smooth data
    Bad = find(isnan(ToPlot)); ToPlot = inpaint_nans(ToPlot);
    ToPlot = smoothn(ToPlot,Settings.SmoothSize); ToPlot(Bad) = NaN;    
    Bad = find(isnan(Wind)); Wind = inpaint_nans(Wind);
    Wind = smoothn(Wind,Settings.SmoothSize); Wind(Bad) = NaN;    
    
    %compute colour and line levels
    ColourLevels = linspace(ColourRange(1),ColourRange(2),16);
    
    %plot data
    ToPlot(ToPlot < min(ColourLevels)) = min(ColourLevels); %to avoid white patches which could be confused with no-data
    contourf(x,y,ToPlot,ColourLevels,'edgecolor','none'); shading flat; hold on

    
    %plot wind
    [c,h] = contour(x,y,Wind,10:10:100,'-','color','k','linewi',2);
    clabel(c,h,'color','k');
    
% %     
% %     if iBasis == 1;
% %       %plot tropause
% %       plot(Lat,Plot.TP(iQuarter,:),'k-','linewi',2)
% %     end

    %tidy
    caxis(ColourRange)
    colormap(cbrew(Settings.ColourTable,16))
    set(gca,'ydir','reverse') %pressure drops with height
    if iBasis == 1; set(gca,'yscale','log'); end
    axis([min(Lat)-5 max(Lat)+5 min(Scale) max(Scale)])
    
% %     %grid lines
% %     if iBasis == 1;
% %       for y = 300:-50:-200; plot([-100,100],[1,1].*y,'--','color',[1,1,1].*0.6); end
% %       for x = -80:20:80;   plot([1,1].*x,[1e-5,500],'--','color',[1,1,1].*0.6); end
% %       plot([0,0],[1e-5,500],'k--','linewi',1)
% %     else
% % 
% %       for y = -400:50:250; plot([-100,100],[1,1].*y,'--','color',[1,1,1].*0.6); end      
% %       for x = -80:20:80;  plot([1,1].*x,[-500,500],'--','color',[1,1,1].*0.6); end
% %       plot([0,0],[-500,500],'k--','linewi',1)
% %       plot([-100,100],[0,0],'k-','linewi',2)
% %     end
    
    
    
    %labelling
    if     iQuarter == 4; set(gca,'yaxislocation','right'); 
    elseif iQuarter ~= 1; set(gca,'yticklabel',[]); 
    end
    if iQuarter == 1 | iQuarter == 4; ylabel(YLabel); end
    
    if iBasis == 1; set(gca,'xaxislocation','top'); end
    xlabel('Latitude [deg]')
    
    if iBasis == 1; text(-39, 330,['(',Letters(k),')'],'fontsize',20); end    
    if iBasis == 2; text(-39,-370,['(',Letters(k),')'],'fontsize',20); end  
    
    if iBasis == 2;
      switch iQuarter;
        case 1; title('DJF','fontsize',20);
        case 2; title('MAM','fontsize',20);
        case 3; title('JJA','fontsize',20);
        case 4; title('SON','fontsize',20);
      end
    end
    
%     if iBasis == 2; ylim([-300 300]); end
          
    drawnow
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% overplot winds
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    
    
  end
end

%colourbar
cb =colorbar('position',[0.92 0.33 0.02 0.33]);
cb.Label.String = 'Gravity Wave Amplitude [K]';
  