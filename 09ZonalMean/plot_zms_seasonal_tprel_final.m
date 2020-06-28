clearvars


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot seasonal zonal mean-height plots
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/06/03
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Settings.InFile = 'zm_final_lonweighted.mat';

Settings.Seasons = [12,1,2;3,4,5;6,7,8;9,10,11];

Settings.SmoothSize = [1,1].*3; %in raw data bins
Settings.InterpFactor = [5,5]; 
    
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prepare data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load file
Data = load(Settings.InFile);

%pull out scales and titles
Lat   = Data.Settings.LatScale;
Time  = Data.Settings.TimeScale;
Scale = Data.Settings.RelPrsScale;

%generate plot data
Plot.GWs  = NaN(4,numel(Lat),numel(Scale));
Plot.Wind = NaN(4,numel(Lat),numel(Scale));
Plot.N    = NaN(4,numel(Lat),numel(Scale));
Plot.k    = NaN(4,numel(Lat),numel(Scale));
Plot.TP   = NaN(4,numel(Lat));

%months, for use in season identification
[~,mm,~] = datevec(Time);

for iSeason = 1:1:4;
  
  %identify days in season
  ThisSeason = find(ismember(mm,Settings.Seasons(iSeason,:)));
  
  %average wrt time and store
  Plot.GWs( iSeason,:,:) = squeeze(nanmean(Data.Results.Rel( 2,ThisSeason,:,:),2));
  Plot.Wind(iSeason,:,:) = squeeze(nanmean(Data.Results.Rel( 1,ThisSeason,:,:),2));  
  Plot.N(   iSeason,:,:) = squeeze(nansum( Data.Results.RelN(2,ThisSeason,:,:),2));  
  Plot.k(   iSeason,:,:) = 1./squeeze(nansum(Data.Results.Rel( 3,ThisSeason,:,:),2));
  Plot.k(isinf(Plot.k))= NaN;

  
  
end; clear iSeason mm Time ThisSeason 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% special handling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for N, log the results and is also a function of smoothsize
Plot.N = log10(Plot.N .* prod(Settings.SmoothSize));
Plot.N(isinf(Plot.N)) = NaN;

%rescale smoothing for overinterpolation
Settings.SmoothSize = Settings.InterpFactor.*Settings.SmoothSize;
if ~isodd(Settings.SmoothSize(1)); Settings.SmoothSize(1) = Settings.SmoothSize(1)+1; end
if ~isodd(Settings.SmoothSize(2)); Settings.SmoothSize(2) = Settings.SmoothSize(2)+1; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf
set(gcf,'color','w')
subplot = @(m,n,p) subtightplot (m, n, p, [0.02, 0.02], 0.08, [0.06 0.09]);

Letters = 'abcdefghjklmnopqrstuvwxyz';
k = 0;

for iVariable=1:1:4;
  
  switch iVariable
    case 1; 
      Data = Plot.GWs;   
      ColourLevels =  0.3:0.04:1.2;
      LineLevels   =  0:0.25:1.5;
      Colours = cbrew('RdYlBu',numel(ColourLevels)+1);
      Name = 'GW Amplitude [K]';
    case 2; 
      Data = Plot.Wind;  
      ColourLevels = -35:2.5:35;
      LineLevels   =  -40:10:40;
      Colours = cbrew('nph_BlueOrange',numel(ColourLevels)+1);
      Name = 'Zonal Wind [m/s]';
    case 3; 
      Data = Plot.k;  
      ColourLevels = 0:20:300;
      LineLevels   = 0:60:1000;
      Colours = cbrew('nph_BlueOrange',numel(ColourLevels)+1);
      Name = 'Wavelength [km]';      
    case 4; Data  = Plot.N;
      Data = Plot.N; 
      ColourLevels = 3:.25:7;
      LineLevels   = 0:0.5:7;
      Colours = cbrew('Blues',numel(ColourLevels)+1);
      Name = '[log_{10}(N)]';     
  end
  
  
  
  for iQuarter = 1:1:4;
    
    k = k+1; pp = subplot(4,4,k);
    
    %pull data
    ToPlot = squeeze(Data(iQuarter,:,:))';
    
    %overinterpolate the data. This is so that pcolor/countourf show
    %data when one vertex is NaN
    x = linspace(  min(Lat),  max(Lat),ceil(range(  Lat)./(mean(diff(  Lat))./Settings.InterpFactor(1))));
    y = linspace(min(Scale),max(Scale),ceil(range(Scale)./(mean(diff(Scale))./Settings.InterpFactor(2))));
    [xi,yi] = meshgrid(x,y);
    ToPlot = interp2(Lat,Scale,ToPlot,xi,yi,'nearest');
    
    %smooth data
    Bad = find(isnan(ToPlot)); ToPlot = inpaint_nans(ToPlot);
    ToPlot = smoothn(ToPlot,Settings.SmoothSize); ToPlot(Bad) = NaN;    
    
    
    %plot data
    ToPlot(ToPlot < min(ColourLevels)) = min(ColourLevels); %to avoid white patches which could be confused with no-data
    contourf(x,y,ToPlot,ColourLevels,'edgecolor','none'); shading flat; hold on

    %line contours
    [c,h] = contour(x,y,ToPlot,LineLevels,'edgecolor',[1,1,1].*0.3); 
    clabel(c,h)
    
    
    %plot tropause
    plot(Lat,zeros(size(Lat)),'k-','linewi',3)

    %tidy
    caxis([min(ColourLevels) max(ColourLevels)])
    colormap(pp,Colours)
    set(gca,'ydir','reverse') %pressure drops with height
%     set(gca,'yscale','log');
    axis([min(Lat)-2.5 max(Lat)+2.5 -250 150])
    
    %grid lines
    for y = 300:-50:-200; plot([-100,100],[1,1].*y,':','color',[1,1,1].*0.8,'linewi',0.25); end
    for x = -80:20:80;    plot([1,1].*x,[1e-5,500],':','color',[1,1,1].*0.8,'linewi',0.25); end
    plot([0,0],[-500,500],'k--','linewi',1)

    
    
    %labelling
%     if     iQuarter == 4; set(gca,'yaxislocation','right'); 
    if iQuarter ~= 1; set(gca,'yticklabel',[]); 
    end
    if iQuarter == 1; ylabel('Pressure [hPa]'); end
    
    if iVariable ~= 4; set(gca,'xticklabel',{});
    else; xlabel('Latitude [deg]'); end
    
    text(-39, 325,['(',Letters(k),')'],'fontsize',20);
    
    
    if iVariable == 1;
      switch iQuarter;
        case 1; title('DJF','fontsize',30);
        case 2; title('MAM','fontsize',30);
        case 3; title('JJA','fontsize',30);
        case 4; title('SON','fontsize',30);
      end
    end
          
    if iQuarter == 1
      switch iVariable
        case 1; y = 0.75;
        case 2; y = 0.50;
        case 3; y = 0.25;
        case 4; y = 0.00;          
      end
      cb =colorbar('position',[0.92 y 0.02 0.2]);
      cb.Label.String = Name;
    end
    
    drawnow
    
  end
end

