clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot time series of GW properties measured by IAGOS
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/Jun/07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Settings.InFile = 'northatlantic_31.mat';%'newfoundland_14_b80.mat';
% Settings.InFile = 'greenland_31_b80.mat'; 

%plots
Settings.Combos = {{'STT_A',4,'N. Flights',[0 40]},      ... 
                   {'STT_A',3,'Amplitude',[0.2,1.8]},      ... 
                   {'STT_k',3,'Wavelength',[50 500]}};
%                    {'T',2,'Temperature',[210 230]}};            %mean temperature
                 %{'U',3,'Median wind speed',[0,60]}
%                    {'STT_A',1,'N. points [x10^4]',[0,2.25]},  ... %number of points used for amplitude

Settings.SmoothSize = 7;
Settings.LongSmooth = 61;

Letters = 'abcdefghhijklmnopqrstuvwxyz'; k =0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prep page
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf
set(gcf,'color','w')
subplot = @(m,n,p) subtightplot (m, n, p, 0.05, 0.05, 0.1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load and prep data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load file
Data = load(Settings.InFile);

Data.Results.STT_k = 1./Data.Results.STT_k;
Data.Results.STT_A(:,1) = Data.Results.STT_A(:,1)./1e4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot data - long time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iCombo=1:1:numel(Settings.Combos);

  %generate subplot panel
  subplot(numel(Settings.Combos),7,[1:5]+(iCombo-1).*7)
  
  %extract time series
  Combo = Settings.Combos{iCombo};
  TS     = Data.Results.(Combo{1});
  TS     = TS(:,Combo{2});
  Bad = find(isnan(TS)); TS= inpaint_nans(TS);
  LTS    = smoothn(TS,[Settings.LongSmooth,1]);
  TS     = smoothn(TS,[Settings.SmoothSize,1]);
  LTS(Bad) = NaN; TS(Bad) = NaN;
  
  %shade years
  for iYear=1994:2:2020;
    y = Combo{4}; yp = [y(1),y(1),y(2),y(2),y(1)];
    patch(datenum([iYear,iYear+1,iYear+1,iYear,iYear],1,1),yp,[1,1,1].*0.9,'edgecolor','none')
    hold on
  end
  
  %plot data
  plot(Data.Settings.TimeScale,TS,'-','color',[1,1,1].*0.6,'linewi',0.5); 
  plot(Data.Settings.TimeScale,LTS,'k-','linewi',2); 
  
  %linear trend
  Good = find(~isnan(Data.Settings.TimeScale + LTS'));
  [p,S] = polyfit(Data.Settings.TimeScale(Good),LTS(Good),1);
  px = polyval(p,Data.Settings.TimeScale);
  plot(Data.Settings.TimeScale,px,'r-')
  
  %tidy
  datetick
  if     iCombo == 1;                     set(gca,'xaxislocation','top');
  elseif iCombo ~= numel(Settings.Combos);set(gca,'xticklabel',{});
  end
  
  ylabel(Combo{3})
  ylim(Combo{4})
  xlim([datenum(1994,8,1),datenum(2019,12,31)])
  box on; grid off

  set(gca,'xtick',datenum([1994:2:2020],7,1),'xticklabel',[1994:2:2020],'fontsize',10)  
  k = k+1; text(datenum(1994,8,5),y(1)+0.92.*(y(2)-y(1)),['(',Letters(k),')'],'fontsize',24)
  drawnow
  
end
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot data - annualised series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iCombo=1:1:numel(Settings.Combos);

  %generate subplot panel
  subplot(numel(Settings.Combos),7,[-1,0]+iCombo.*7)

  
  
  %extract time series
  Combo = Settings.Combos{iCombo};
  TS = Data.Results.(Combo{1});
  TS = TS(:,Combo{2});
  TS = smoothn(TS,Settings.SmoothSize);
  
  %plot each year
  [yy,~,~] = datevec(Data.Settings.TimeScale);
  years = unique(yy); years = years(~isnan(years));
  dd = date2doy(Data.Settings.TimeScale);
  for iYear=1:1:numel(years)
    ThisYear = find(yy == years(iYear));
    plot(dd(ThisYear),TS(ThisYear),'-','color',[1,1,1].*0.6,'linewi',0.5)
    hold on
  end
  
  %overplot all years mean
  All = NaN(365,1);
  for iDay=1:1:365;
    All(iDay) = nanmean(TS(dd == iDay));
  end
  
  plot(1:1:365,All,'k-','linewi',3)
  
%   %plot data
%   plot(Data.Settings.TimeScale,TS)
%   
  %tidy
  set(gca,'yaxislocation','right')
  if     iCombo == 1;                     set(gca,'xaxislocation','top');
  elseif iCombo ~= numel(Settings.Combos);set(gca,'xticklabel',{});
  end  
  ylabel(Combo{3}); set(gca,'yaxislocation','right')
  datetick('x','m'); grid off
  set(gca,'xtick',datenum(0,1:1:12,1)+15)  
  ylim(Combo{4})
  
  y = Combo{4};
  k = k+1; text(5,y(1)+0.92.*(y(2)-y(1)),['(',Letters(k),')'],'fontsize',24)
  drawnow
  
end
  
