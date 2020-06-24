clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot IAGOS postprocessed metadata - time distribution
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/05/14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%file handling
Settings.InFile = 'metadata_all_sgolay900.mat';

%height regions (absolute pressure)
Settings.Region.Bot = [1000,225];
Settings.Region.Mid = [225,205];
Settings.Region.Top = [205,0];
Settings.Region.All = [1000,0];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load data and split by region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load data
Data = load(Settings.InFile);

%define regions
idx.Bot = find(Data.Settings.Grid.Prs >= min(Settings.Region.Bot) & Data.Settings.Grid.Prs <= max(Settings.Region.Bot));
idx.Mid = find(Data.Settings.Grid.Prs >= min(Settings.Region.Mid) & Data.Settings.Grid.Prs <= max(Settings.Region.Mid));
idx.Top = find(Data.Settings.Grid.Prs >= min(Settings.Region.Top) & Data.Settings.Grid.Prs <= max(Settings.Region.Top));
idx.All = find(Data.Settings.Grid.Prs >= min(Settings.Region.All) & Data.Settings.Grid.Prs <= max(Settings.Region.All));

%hence, make three calendars
Calendar.Bot = nansum(Data.Results.Time.Prs(:,:,idx.Bot),3);
Calendar.Mid = nansum(Data.Results.Time.Prs(:,:,idx.Mid),3);
Calendar.Top = nansum(Data.Results.Time.Prs(:,:,idx.Top),3);
Calendar.All = nansum(Data.Results.Time.Prs(:,:,idx.All),3);
sum(Calendar.All(:))
% stop
Years  = Data.Settings.Grid.Years;
Months = Data.Settings.Grid.Months;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define colours
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Range = [0, prctile([Calendar.Bot(:);Calendar.Mid(:);Calendar.Top(:)],85)];
Range = ceil(Range./10000).*10000;
Settings.NColours = ceil(max(Range)./10000).*2;

Levels = linspace(Range(1),Range(2),Settings.NColours);

Colours = flipud(cbrew('RdYlGn',Settings.NColours+1));
Colours = Colours(1:end-1,:); %The really dark green is hard to read text over


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot data. do this by hand with patches, for more control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%prepare figure
clf
set(gcf,'color','w')
subplot = @(m,n,p) subtightplot (m, n, p, 0.028, 0.1, [0.03,0.08]);

%plot data
for iRegion=1:1:4;
  
  %create panel
  subplot(1,4,iRegion)
  axis([1 13 1994 2020])  
  
  %get data
  switch iRegion
    case 1; Cal = Calendar.Bot; title('(a) Below 225 hPa');
    case 2; Cal = Calendar.Mid; title('(b) 225 hPa - 205 hPa');
    case 3; Cal = Calendar.Top; title('(c) Above 205 hPa');      
    case 4; Cal = Calendar.All; title('(d) All Data');       
  end
  
  
  %plot data  
  for iYear=1:1:numel(Years)-1
    for iMonth = 1:1:numel(Months)
      
      %define the patch limits
      x = [0,1,1,0,0]+iMonth;
      y = [0,0,1,1,0]+Years(iYear);
      
      %work out the colour
      N = Cal(iYear,iMonth);
      idx = closest(Levels,N);
      
      
      %plot the patch
      if N == 0; patch(x,y,'w','edgecolor',[1,1,1].*0.6);
      else
        
        if  iRegion ~=4; Colour = Colours(idx,:);
        else;            Colour = [1,1,1].*0.7;
        end
        
        patch(x,y,Colour,'edgecolor',[1,1,1].*0.6)
        
        %label with number of points
        PrintN = round(N./1000);%round(log10(N));
        text(iMonth+0.5,Years(iYear)+0.5,num2str(PrintN), ...
             'verticalalignment','middle','horizontalalignment','center', ...
             'color','k','fontweight','normal','fontsize',12)
      end
    end
    
    %sum for each row
    Sigma = round(sum(Cal(iYear,:))./1000);
    text(13.2,Years(iYear)+0.5,num2str(Sigma), ...
         'verticalalignment','middle','horizontalalignment','left', ...
         'color','k','fontweight','normal','clipping','off','fontsize',14)
  end
  
  %sum for each column
  for iMonth=1:1:12;
    Sigma = round(sum(Cal(:,iMonth))./1000);
    if mod(iMonth,2) == 0; Style = '\it'; else Style = '\rm'; end
    
    text(iMonth+0.43,2020.5,[Style,'{',num2str(Sigma),'}'], ...
         'verticalalignment','middle','horizontalalignment','center', ...
         'color','k','fontweight','normal','clipping','off','fontsize',10)
% %     if iMonth < 12;
% %       text(iMonth+1,2021.5,'¦', ...
% %         'verticalalignment','middle','horizontalalignment','center', ...
% %         'color','k','fontweight','bold','clipping','off','fontsize',12)
% %     end
  end
  
  

  

  
  %tidy
  if iRegion== 1;   set(gca,'ytick',Years+0.5,'yticklabel',Years,'xaxislocation','top')
  else                 set(gca,'ytick',[],'yticklabel',{},'xaxislocation','top')
  end
  set(gca,'ydir','reverse','fontweight','bold')
  set(gca,'xtick',1.5:1:12.5, 'xticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D'})

  drawnow
      
end



%% done! colourbar
colormap(Colours)
cb = colorbar('position',[0.952,0.25,0.01,0.5],'fontweight','normal');
caxis(Range./1000)
cb.Label.String = 'Thousands of observations';


