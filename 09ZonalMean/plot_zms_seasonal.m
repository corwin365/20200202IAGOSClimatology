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

Settings.InFile = 'zm_longrange.mat';

Settings.Variable = 3; %in order of input file
Settings.Mode = 1; %1 for results, 2 for number of points

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
Var     = Data.Settings.Vars{Settings.Variable};

%generate plot data
Plot.Abs = NaN(4,numel(Lat),numel(Prs.Abs));
Plot.Rel = NaN(4,numel(Lat),numel(Prs.Rel));
Plot.TP  = NaN(4,numel(Lat));

%months, for use in season identification
[~,mm,~] = datevec(Time);

for iSeason = 1:1:4;
  
  %identify days in season
  ThisSeason = find(ismember(mm,Settings.Seasons(iSeason,:)));
  
  %average wrt time and store
  if Settings.Mode == 1;
    Plot.Abs(iSeason,:,:) = squeeze(nanmean(Data.Results.Abs(Settings.Variable,ThisSeason,:,:),2));
    Plot.Rel(iSeason,:,:) = squeeze(nanmean(Data.Results.Rel(Settings.Variable,ThisSeason,:,:),2));
  elseif Settings.Mode == 2;
    Plot.Abs(iSeason,:,:) = squeeze(nanmean(Data.Results.AbsN(Settings.Variable,ThisSeason,:,:),2));
    Plot.Rel(iSeason,:,:) = squeeze(nanmean(Data.Results.RelN(Settings.Variable,ThisSeason,:,:),2));
  end

  %tropopause height is done separately either way
  Plot.TP(iSeason,:) = nanmean(Data.Results.TP(ThisSeason,:));
  
  
  
end; clear iSeason mm Time ThisSeason 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% special handling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%k is upside down
if Settings.Mode == 1 && strcmp(Var,'STT_k'); Plot.Abs = 1./Plot.Abs; Plot.Rel = 1./Plot.Rel; end

%colour table
if Settings.Mode == 2; Settings.ColourTable = 'Reds';
else
  switch Var
    case 'T';     Settings.ColourTable = 'RdYlBu';
    case 'U';     Settings.ColourTable = 'RdBu';
    case 'STT_A'; Settings.ColourTable = 'nph_BlueOrange';
    case 'STT_k'; Settings.ColourTable = 'PRGn';
  end
end

%for N, log the results
if Settings.Mode == 2; 
  Plot.Abs = log10(Plot.Abs);      Plot.Rel = log10(Plot.Rel);  
  Plot.Abs(isinf(Plot.Abs)) = NaN; Plot.Rel(isinf(Plot.Rel)) = NaN;
end

%for N, variable name is slightly different
if Settings.Mode == 2; Var = 'log(Number of points)'; end

%colour range
if     Settings.Mode == 1; 
  if strcmp(Var,'U') == 1; ColourRange = [-1,1].*prctile(abs([Plot.Abs(:);Plot.Rel(:)]),[95]);
  else;                    ColourRange = prctile([Plot.Abs(:);Plot.Rel(:)],[5,95]);
  end
elseif Settings.Mode == 2; ColourRange = [0,max([Plot.Abs(:);Plot.Rel(:)])];
end

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
    case 1; Data = Plot.Abs; Scale = Prs.Abs; YLabel = 'Pressure [hPa]';
    case 2; Data = Plot.Rel; Scale = Prs.Rel; YLabel = '\Delta Pressure [hPa]';
  end
  
  for iQuarter = 1:1:4;
    
    k = k+1; subplot(2,4,k);
    
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
    
    
    %compute colour and line levels
    ColourLevels = linspace(ColourRange(1),ColourRange(2),16);
    
    %plot data
    ToPlot(ToPlot < min(ColourLevels)) = min(ColourLevels); %to avoid white patches which could be confused with no-data
    contourf(x,y,ToPlot,ColourLevels,'edgecolor','none'); shading flat; hold on

    
    if iBasis == 1;
      %plot tropause
      plot(Lat,Plot.TP(iQuarter,:),'k-','linewi',2)
    end

    %tidy
    caxis(ColourRange)
    colormap(cbrew(Settings.ColourTable,16))
    set(gca,'ydir','reverse') %pressure drops with height
    if iBasis == 1; set(gca,'yscale','log'); end
    axis([min(Lat)-5 max(Lat)+5 min(Scale) max(Scale)])
    
    %grid lines

    
    if iBasis == 1;
      for y = 300:-50:-200; plot([-100,100],[1,1].*y,'--','color',[1,1,1].*0.6); end
      for x = -80:20:80;   plot([1,1].*x,[1e-5,500],'--','color',[1,1,1].*0.6); end
      plot([0,0],[1e-5,500],'k--','linewi',1)
    else

      for y = -400:50:250; plot([-100,100],[1,1].*y,'--','color',[1,1,1].*0.6); end      
      for x = -80:20:80;  plot([1,1].*x,[-500,500],'--','color',[1,1,1].*0.6); end
      plot([0,0],[-500,500],'k--','linewi',1)
      plot([-100,100],[0,0],'k-','linewi',2)
    end
    
    
    
    %labelling
    if     iQuarter == 4; set(gca,'yaxislocation','right'); 
    elseif iQuarter ~= 1; set(gca,'yticklabel',[]); 
    end
    if iQuarter == 1 | iQuarter == 4; ylabel(YLabel); end
    
    if iBasis == 1; set(gca,'xaxislocation','top'); end
    xlabel('Latitude [deg]')
    
    if iBasis == 1; text(-39, 178,['(',Letters(k),')'],'fontsize',20); end    
    if iBasis == 2; text(-39,-370,['(',Letters(k),')'],'fontsize',20); end  
    
    if iBasis == 2;
      switch iQuarter;
        case 1; title('DJF','fontsize',20);
        case 2; title('MAM','fontsize',20);
        case 3; title('JJA','fontsize',20);
        case 4; title('SON','fontsize',20);
      end
    end
    
    if iBasis == 2; ylim([-200 200]); end
          
    drawnow
    
  end
end

%colourbar
cb =colorbar('position',[0.92 0.33 0.02 0.33]);
cb.Label.String = Var;
  