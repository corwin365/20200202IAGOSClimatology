clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%plot individual time series for each identified peak
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %% part A: draw map
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % clf
% % set(gcf,'color','w')
% % 
% % %generate projection
% % %%%%%%%%%%%%%%%%%%%%%%
% % 
% % m_proj('lambert','lat',[25,80],'lon',[-130 150])
% % 
% % %plot land
% % %%%%%%%%%%%%%%%%%%%%%%%
% % m_coast('patch',[178,255,102]./255,'edgecolor',[1,1,1].*0.5); %coasts at the bottom
% % 
% % %plot locations
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%

Letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
Lons = [-115, -105, -83,  -83,-90,-73,-69, -60, -47, -18, -8,  9, 12, 25, 35, 49, 60, 67, 69, 100, 90, 120, 138, 117];
Lats = [  45,   52,  32, 46.5, 59, 54, 68,  48,  64,  65, 59, 56, 51, 35, 49, 35, 62, 45, 38,  65, 52,  55,  48,  38];
Oro =  [1,0,0.5,0,0,0.5,1,0,1,1,0.5,0,0.5,0,0,1,1,1,1,1,1,1,1,0];
Names = {'Rockies','Canadian Plains','Georgia, US','Great Lakes','Hudson Bay','Quebec','Baffin Island','Newfoundland','Greenland','Iceland','Scotland','South Scandinavia','European Plain','East Med','Ukraine','Iran','Urals','South Kazakhstan','Afghanistan','Siberia','Altai','CMR Corner','Sikhote Alin','East China'};

% % Order = [8;16;13;2;4;19;14;12;23;24;22;1;18;21;15;11;3;10;5;9;20;17;6;7]; %determined below, and hardcoded
% % Lons = Lons(Order);
% % Lats = Lats(Order);
% % Oro = Oro(Order);
% % Names = Names(Order);

% % for iMarker=1:1:numel(Lats)
% %   m_text(Lons(iMarker),Lats(iMarker),Letters(iMarker),'fontsize',18,'fontweight','bold')
% % end; clear iMarker
% % 
% % %gridlines and done
% % %%%%%%%%%%%%%%%%%%%%%%%
% % m_grid('fontsize',18,'linestyle','--','color',[1,1,1].*.4);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate the time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%what data do we want?
Range = [10,5]; %window size for average
Generate = 0; %only generate new data if needed, it is *very* slow


if Generate == 1;
  
  for iSeries = 1:1:numel(Lats);
    disp(['Processing box ',Letters(iSeries)])
    
    
    LatRange = 0.5.*Range(2)*[-1,1] + Lats(iSeries);
    LonRange = 0.5.*Range(1)*[-1,1] + Lons(iSeries);
    
    
    func_generate_time_series(Letters(iSeries), LonRange, LatRange, 'STT_A', 'nanmedian', 31)
  end
end
clear LatRange LonRange iSeries Generate Range

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% annualise the time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Series = NaN(numel(Lats),366,3); %3 is [18,50,82]th percentile

for iSeries = 1:1:numel(Lats)
  
  File = ['data/',Letters(iSeries),'.mat'];
  if ~exist(File,'file'); continue; end
  Data = load(File);
  dd = date2doy(Data.Settings.TimeScale);
  for iDay=1:1:366;
    OnThisDay = find(dd == iDay);
    DayData = Data.Results.('STT_A');
    DayData = DayData(OnThisDay);
    Good = find(~isnan(DayData));
    if numel(Good) ~= 0; Series(iSeries,iDay,:) = prctile(DayData(Good),[18,50,82]); end
  end
end
clear iSeries File Data dd iDay OnThisDay DayData Good

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %% sort by similarity to Newfoundland (arbitrary choice, cleanest by eye)
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % %metric one: correlation
% % Metric1 = NaN(numel(Lats),1);
% % for iM=1:1:numel(Lats)
% %   s1 = Series(8,:,2);
% %   s2 = Series(iM,:,2);
% %   Good = find(~isnan(s1+s2));
% %   cc = corrcoef(s1(Good),s2(Good));
% %   Metric1(iM) = cc(2);
% % end
% % 
% % %metric two: rms difference
% % Metric2 = NaN(numel(Lats),1);
% % for iM=1:1:numel(Lats)
% %   s1 = Series(8,:,2);
% %   s2 = Series(iM,:,2);
% %   Good = find(~isnan(s1+s2));
% %   Metric2(iM) = sqrt(sum((s1(Good)-s2(Good)).^2));
% % end
% % 
% % 
% % %combine
% % Sigma = Metric1 + 1-(Metric2./nanmax(Metric2));
% % [~,Order] = sort(Sigma,'descend');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot them as small multiples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


SmoothSize = 13;

clf
set(gcf,'color','w')
subplot = @(m,n,p) subtightplot (m, n, p, [0.012, 0.045], [0.05 0.05], [0.25 0.2]);

%assign plotting order
Order = 1:1:numel(Lats);

k = 0;
for jSeries = 1:1:numel(Order);
  iSeries = Order(jSeries);
  
  %get the data
  ThisSeries = squeeze(Series(iSeries,:,:));
  ThisSeries = smoothn(ThisSeries,[SmoothSize,1]);
  
  if nansum(ThisSeries(:)) == 0; continue; end
  
  %find the series median
  SeriesMedian = nanmedian(ThisSeries(:,2)); 
  
  %and a range to plot over
  YRange = [nanmin(ThisSeries(:,2)),nanmax(ThisSeries(:,2))]+[-1,1].*0.1;
  
  %generate subplot
  k = k+1;
  subplot(ceil(numel(Lats)./3),3,k)
  axis([0 366 YRange]); box on, grid off;
  hold on
  
  %series labels
  text(-50,mean(YRange),Letters(jSeries),'fontsize',24,'fontweight','bold',...
      'verticalalignment','middle','horizontalalignment','right');

  
  %axis tick locations
  if     k < 4;               set(gca,'xaxislocation','top');
  elseif k > numel(Series)-3; set(gca,'xaxislocation','bottom');
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
  if  k < 4 | k > numel(Series)-3; set(gca,'xtick',Days,'xticklabels',Months);
  else;                            set(gca,'xtick',Days,'xticklabels',{});
  end
  set(gca,'tickdir','out')
  
  %plot a line under the data just so we don't get discontinuities at section edges
  plot(1:1:366,ThisSeries(:,2),'b-','linewi',1) 
  
  
  %plot above-average data
  x = 1:1:366;
  Above = find(ThisSeries(:,2) >= SeriesMedian);
  ToPlot.x = [1:1:366,366:-1:1,1];
  ToPlot.y = ones(366,1).*SeriesMedian;
  ToPlot.y(Above) = ThisSeries(Above,2);
  ToPlot.y = [ToPlot.y;ones(366,1).*SeriesMedian;ToPlot.y(1)]';
  patch(ToPlot.x,ToPlot.y,[204,229,255]./255,'edgecolor','none') 
  
  %plot below-average data
  x = 1:1:366;
  Below = find(ThisSeries(:,2) <=SeriesMedian);
  ToPlot.x = [1:1:366,366:-1:1,1];
  ToPlot.y = ones(366,1).*SeriesMedian;
  ToPlot.y(Below) = ThisSeries(Below,2);
  ToPlot.y = [ToPlot.y;ones(366,1).*SeriesMedian;ToPlot.y(1)]';
  patch(ToPlot.x,ToPlot.y,[255,204,204]./255,'edgecolor','none') 
  
  %plot data limit
  x2 = x; x2(Below) = NaN;
  y2 = ThisSeries(:,2); y2(Below) = NaN;
  plot(x2,y2,'b-','linewi',1) 
  
  x2 = x; x2(Above) = NaN;
  y2 = ThisSeries(:,2); y2(Above) = NaN;
  plot(x2,y2,'r-','linewi',1)   

  
  %median line
  plot([0,366],[1,1].*SeriesMedian,'color',[1,1,1].*0.3)

  %mountains nearby?
  if Oro(jSeries) == 1;
    text(330,YRange(1)+0.85*range(YRange),'\it\Lambda','fontweight','bold','fontsize',24);
  elseif Oro(jSeries) == 0.5
    text(330,YRange(1)+0.87*range(YRange),'\it\Lambda','fontsize',16);
  end
  
  %names
  if nanmean(ThisSeries(1:90,2)) > SeriesMedian
    text(2,YRange(1)+0.1*range(YRange),Names{jSeries},'fontweight','normal','fontsize',12);
  else
    text(2,YRange(1)+0.9*range(YRange),Names{jSeries},'fontweight','normal','fontsize',12);
  end
  
  %done
  drawnow
  
end


% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %% additional plot: absolute range
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % % % figure
% % % clf
% % % axis([0.5 24.5 0.05 1.6])
% % % hold on; box on; grid off;
% % % 
% % % 
% % % %horizontal lines
% % % for iY=0.2:0.2:1.4; plot([0.5,24.5],[1,1].*iY,'--','color',[1,1,1].*0.8); end
% % % set(gca,'ytick',0.2:0.2:1.4)
% % % 
% % % for iSeries = 1:1:24;
% % % 
% % %   jSeries = Order(iSeries);
% % %   
% % %   %background line
% % %   plot([1,1].*iSeries,[0 2.2],'-','color',[1,1,1].*0.8)
% % %   
% % %   %get data
% % %   Data = smoothn(squeeze(Series(jSeries,:,2)),[SmoothSize,1]);
% % %   Percentiles = prctile(Data(~isnan(Data)),[0,2.5,18,50,82,97.5,100]);
% % % 
% % %   %plot stats
% % %   %%%%%%%%%%%%%%
% % % %   %min and max
% % % %   plot(iSeries,Percentiles(1),'r.','markersize',12,'linewi',2)
% % % %   plot(iSeries,Percentiles(7),'b.','markersize',12,'linewi',2)
% % %   
% % %   %2 st dev
% % %   plot(0.2.*[-1,1] + iSeries,[1,1].*Percentiles(2),'r-','linewi',2)  
% % %   plot(0.2.*[-1,1] + iSeries,[1,1].*Percentiles(6),'b-','linewi',2)    
% % %   
% % %   %-1 st dev to +1 stdev
% % %   patch(iSeries+[-1,1,1,-1,-1].*0.2,Percentiles([3,3,4,4,3]),[255,204,204]./255,'edgecolor','r')
% % %   patch(iSeries+[-1,1,1,-1,-1].*0.2,Percentiles([4,4,5,5,4]),[204,229,255]./255,'edgecolor','b')
% % %   
% % %   %median
% % %   plot(0.2.*[-1,1] + iSeries,[1,1].*Percentiles(4),'k-','linewi',2)  
% % %   
% % %   %vertical joins
% % %   plot([1,1].*iSeries,Percentiles([2,3]),'r-')
% % %   plot([1,1].*iSeries,Percentiles([5,6]),'b-')  
% % %   
% % %   %mountains nearby?
% % %   if Oro(iSeries) == 1;
% % %     text(iSeries,0.13,'\it\Lambda','fontweight','bold','fontsize',14,'horizontalalignment','center');
% % %   elseif Oro(iSeries) == 0.5
% % %     text(iSeries,0.13,'\it\Lambda','fontsize',8,'horizontalalignment','center');
% % %   end
% % % end
% % % 
% % % set(gca,'xtick',1:24,'xticklabel',{'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X'})
% % % ylabel('Amplitude [K]')
% % % set(gca,'yaxislocation','right')
% % % 
