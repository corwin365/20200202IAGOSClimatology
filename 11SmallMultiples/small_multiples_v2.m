clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate individual time series for each identified peak
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%time and pressure range
TimeScale = datenum(1994,8,1):1:datenum(2019,12,31);
% TimeScale = datenum(2000,1,1):1:datenum(2005,12,31);
PrsRange = [10000,0];
TimeWindow = 31; %days

%what data do we want?
Range = 500; %km
Generate = 0; %only generate new data if needed, it is *very* slow

%baseline time series
Baseline = 1;

%plot smoothing
SmoothSize = 7;

%individual series to generate
Names = {}; Lons = []; Lats = []; Oro = [];
Names{end+1} = 'Rockies';           Lons(end+1) = -110; Lats(end+1) = 40; 
Names{end+1} = 'Iceland';           Lons(end+1) =  -18; Lats(end+1) = 65; 
Names{end+1} = 'Newfoundland';      Lons(end+1) =  -60; Lats(end+1) = 48; 
Names{end+1} = 'Greenland';         Lons(end+1) =  -47; Lats(end+1) = 64;
Names{end+1} = 'UK';                Lons(end+1) =   -2; Lats(end+1) = 54; 
Names{end+1} = 'Iran';              Lons(end+1) =   49; Lats(end+1) = 35;
Names{end+1} = 'Altai';             Lons(end+1) =   90; Lats(end+1) = 52;
Names{end+1} = 'Sikhote Alin';      Lons(end+1) =  138; Lats(end+1) = 48;
Names{end+1} = 'Urals';             Lons(end+1) =   65; Lats(end+1) = 62;
Names{end+1} = 'Alps/Balkans';      Lons(end+1) =   15; Lats(end+1) = 46;
Names{end+1} = 'Canadian Plains';   Lons(end+1) = -110; Lats(end+1) = 61;
Names{end+1} = 'Siberia';           Lons(end+1) =  100; Lats(end+1) = 65;
Names{end+1} = 'Great Lakes';       Lons(end+1) =  -83; Lats(end+1) = 46.5;
Names{end+1} = 'North Atlantic';    Lons(end+1) =  -30; Lats(end+1) = 52;
Names{end+1} = 'Azores';            Lons(end+1) =  -25; Lats(end+1) = 37;
Names{end+1} = 'CMR Border';        Lons(end+1) =  120; Lats(end+1) = 55;
Names{end+1} = 'West Russia';       Lons(end+1) =   40; Lats(end+1) = 55; 
Names{end+1} = 'Afghanistan';       Lons(end+1) =   69; Lats(end+1) = 38;

% % % South Scandinavia  9  56
% % % East Med  25  35
% % % South Kazakhstan  67  45
% % % Afghanistan  69  38


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate the time series (if needed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Generate == 1;
  
  for iSeries = 1:1:numel(Names);
    disp(['Processing box over ',Names{iSeries}])
    
    %generate the time series
    func_generate_composite_time_series(urlencode(Names{iSeries}),Lons(iSeries),Lats(iSeries),Range, PrsRange,TimeScale, 'STT_A', 'nanmean', TimeWindow)

    %open the file we just created and store the location metadata, to allow us to decouple the processing
    File = load(['data/',urlencode(Names{iSeries}),'.mat']);
    File.Lon = Lons(iSeries); 
    File.Lat = Lats(iSeries); 
    File.Oro = Oro(iSeries);
    pause(0.1) %read/write time too fast otherwise
    save(['data/',urlencode(Names{iSeries}),'.mat'],'File')
  end
end
clear LatRange LonRange iSeries Generate Range File Lars Lons Oro PrsRange TimeScale TimeWindow


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% annualise the time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Series = NaN(numel(Names),366,3); %3 is [18,50,82]th percentile
Info   = NaN(numel(Names),3); %metadata

for iSeries = 1:1:numel(Names)
  
  File = ['data/',urlencode(Names{iSeries}),'.mat'];
  if ~exist(File,'file'); continue; end
  Data = load(File); Data = Data.File;
  Info(iSeries,:) = [Data.Lon,Data.Lat,Data.Oro];
  
  for iDay=1:1:366;
    DayData = Data.Results.('STT_A');
    DayData = DayData(iDay);
    Good = find(~isnan(DayData));
    if numel(Good) ~= 0; Series(iSeries,iDay,:) = prctile(DayData(Good),[18,50,82]); end
  end
end
clear iSeries File Data dd iDay DayData Good

Lons = Info(:,1); Lats = Info(:,2); Oro = Info(:,3); clear Info

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% choose a time series, then reorder the data in order of similarity
% choose using an index number 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp(['Using ',Names{Baseline},' as baseline and reordering data accordingly']);

%metric one: correlation
Metric1 = NaN(numel(Names),1);
for iM=1:1:numel(Names)
  s1 = Series(Baseline,:,2);
  s2 = Series(iM,:,2);
  Good = find(~isnan(s1+s2));
  cc = corrcoef(s1(Good),s2(Good));
  Metric1(iM) = cc(2);
end

% % %metric two: rms difference
% % Metric2 = NaN(numel(Names),1);
% % for iM=1:1:numel(Names)
% %   s1 = Series(Baseline,:,2);
% %   s2 = Series(iM,:,2);
% %   Good = find(~isnan(s1+s2));
% %   Metric2(iM) = sqrt(sum((s1(Good)-s2(Good)).^2));
% % end


%combine
Sigma = Metric1;% + 0.5.*(1-(Metric2./nanmax(Metric2)));
[~,Order] = sort(Sigma,'descend');

clear Baseline Metric1 iM s1 s2 Good cc Metric2 Sigma

Lats = Lats(Order);
Lons = Lons(Order);
Series = Series(Order,:,:);
Names = Names(Order);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot the individual time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';


clf
set(gcf,'color','w')
subplot = @(m,n,p) subtightplot (m, n, p, [0.012, 0.045], [0.05 0.05], [0.25 0.2]);

for iSeries = 1:1:numel(Names);

  
  %get the data
  ThisSeries = squeeze(Series(iSeries,:,2));
  ThisSeries = smoothn2(ThisSeries',[SmoothSize,1]);
  if nansum(ThisSeries(:)) == 0; continue; end
  
  %find the series median
  SeriesMedian = nanmedian(ThisSeries); 
  
  %and a range to plot over
  YRange = [nanmin(ThisSeries),nanmax(ThisSeries)]+[-1,1].*0.1;
  
  %generate subplot
  subplot(ceil(numel(Lats)./3),3,iSeries)
  axis([0 366 YRange]); box on, grid off;
  hold on

  %series labels
  text(-50,mean(YRange),Letters(iSeries),'fontsize',24,'fontweight','bold',...
      'verticalalignment','middle','horizontalalignment','right');


  %axis tick locations
  if     iSeries < 4;               set(gca,'xaxislocation','top');
  elseif iSeries > numel(Series)-3; set(gca,'xaxislocation','bottom');
  end
  
  %quarter shading
  for iMonth=0:6:12;
    patch([datenum(0,iMonth,1),datenum(0,iMonth+3,-1),datenum(0,iMonth+3,-1),datenum(0,iMonth,1),datenum(0,iMonth,1)], ...
          [YRange(1) YRange(1) YRange(2) YRange(2) YRange(1)], ...
           [1,1,1].*0.9,'edgecolor','none');
  end
  plot([0,366],[1,1].*YRange(1),'k-'); plot([0,366],[1,1].*YRange(2),'k-')
  plot([366,366],YRange,'k-')

  %axis tick values
  Days = date2doy(datenum(0,1:1:12,15)); Months = {'J','F','M','A','M','J','J','A','S','O','N','D'};
  if  iSeries < 4 | iSeries> numel(Series)-3; set(gca,'xtick',Days,'xticklabels',Months);
  else;                                       set(gca,'xtick',Days,'xticklabels',{});
  end
  set(gca,'tickdir','out')

  %plot a line under the data just so we don't get discontinuities at section edges
  plot(1:1:366,ThisSeries,'b-','linewi',1) 
  
  %plot above-average data
  x = 1:1:366;
  Above = find(ThisSeries(:) >= SeriesMedian);
  ToPlot.x = [1:1:366,366:-1:1,1];
  ToPlot.y = ones(366,1).*SeriesMedian;
  ToPlot.y(Above) = ThisSeries(Above);
  ToPlot.y = [ToPlot.y;ones(366,1).*SeriesMedian;ToPlot.y(1)]';
  patch(ToPlot.x,ToPlot.y,[204,229,255]./255,'edgecolor','none') 
  
  %plot below-average data
  x = 1:1:366;
  Below = find(ThisSeries(:) <=SeriesMedian);
  ToPlot.x = [1:1:366,366:-1:1,1];
  ToPlot.y = ones(366,1).*SeriesMedian;
  ToPlot.y(Below) = ThisSeries(Below);
  ToPlot.y = [ToPlot.y;ones(366,1).*SeriesMedian;ToPlot.y(1)]';
  patch(ToPlot.x,ToPlot.y,[255,204,204]./255,'edgecolor','none') 
  
  %plot data limit
  x2 = x; x2(Below) = NaN;
  y2 = ThisSeries(:); y2(Below) = NaN;
  plot(x2,y2,'b-','linewi',1) 
  
  x2 = x; x2(Above) = NaN;
  y2 = ThisSeries(:); y2(Above) = NaN;
  plot(x2,y2,'r-','linewi',1)   

  
  %median line
  plot([0,366],[1,1].*SeriesMedian,'color',[1,1,1].*0.3)

 
    
  %names
  if nanmean(ThisSeries(1:90)) > SeriesMedian
    text(2,YRange(1)+0.1*range(YRange),Names{iSeries},'fontweight','normal','fontsize',12);
  else
    text(2,YRange(1)+0.9*range(YRange),Names{iSeries},'fontweight','normal','fontsize',12);
  end
  
    

  %done
  drawnow
  
  %tidy and loop
  clearvars -except Lats Lons Oro Names Series Letters SmoothSize subplot
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% additional plot: absolute range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
clf
axis([0.5 18.5 0.05 1.6])
hold on; box on; grid off;


%horizontal lines
for iY=0.2:0.2:1.4; plot([0.5,24.5],[1,1].*iY,'--','color',[1,1,1].*0.8); end
set(gca,'ytick',0.2:0.2:1.4)

for iSeries = 1:1:18;


  %background line
  plot([1,1].*iSeries,[0 2.2],'-','color',[1,1,1].*0.8)
  
  %get data
  Data = smoothn(squeeze(Series(iSeries,:,2)),[SmoothSize,1]);
  Percentiles = prctile(Data(~isnan(Data)),[0,2.5,18,50,82,97.5,100]);

  %plot stats
  %%%%%%%%%%%%%%
%   %min and max
%   plot(iSeries,Percentiles(1),'r.','markersize',12,'linewi',2)
%   plot(iSeries,Percentiles(7),'b.','markersize',12,'linewi',2)
  
  %2 st dev
  plot(0.2.*[-1,1] + iSeries,[1,1].*Percentiles(2),'r-','linewi',2)  
  plot(0.2.*[-1,1] + iSeries,[1,1].*Percentiles(6),'b-','linewi',2)    
  
  %-1 st dev to +1 stdev
  patch(iSeries+[-1,1,1,-1,-1].*0.2,Percentiles([3,3,4,4,3]),[255,204,204]./255,'edgecolor','r')
  patch(iSeries+[-1,1,1,-1,-1].*0.2,Percentiles([4,4,5,5,4]),[204,229,255]./255,'edgecolor','b')
  
  %median
  plot(0.2.*[-1,1] + iSeries,[1,1].*Percentiles(4),'k-','linewi',2)  
  
  %vertical joins
  plot([1,1].*iSeries,Percentiles([2,3]),'r-')
  plot([1,1].*iSeries,Percentiles([5,6]),'b-')  
  
end

set(gca,'xtick',1:18,'xticklabel',{'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R'})
ylabel('Amplitude [K]')
set(gca,'yaxislocation','right')

