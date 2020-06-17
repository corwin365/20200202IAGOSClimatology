clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% map indicating locations discussed in text
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% part A: draw map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf
set(gcf,'color','w')

%generate projection
%%%%%%%%%%%%%%%%%%%%%%

m_proj('stereographic','lat',90,'lon',0,'radius',65)
hold on
% [-130 150]

%plot topography
%%%%%%%%%%%%%%%%%%%%%%%

%load image map
[Topo,Map] = topo_etc([-180 179.999],[25,89.9],0,0,0,1);
Map.LonScale(end+1) = 185;
Map.Map(:,end+1,:) = Map.Map(:,1,:);

% % Map.LonScale(Map.LonScale <=0) = Map.LonScale(Map.LonScale <= 0)+360;
% % [~,idx] = sort(Map.LonScale); Map.LonScale = Map.LonScale(idx);
% % Map.Map = Map.Map(:,idx,:);

m_image(Map.LonScale,Map.LatScale,Map.Map);

%plot locations
%%%%%%%%%%%%%%%%%%%%%%%%%%%

Letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
Lons = [-115, -105, -83,  -83,-90,-73,-69, -60, -47, -18, -8,  9, 12, 25, 35, 49, 60, 67, 69, 100, 90, 120, 138, 117];
Lats = [  45,   52,  32, 46.5, 59, 54, 68,  48,  64,  65, 59, 56, 51, 35, 49, 35, 62, 45, 38,  65, 52,  55,  48,  38];

%order is computed in the small_multiples routine, and hardcode here for consistency
Order = [8;16;13;2;4;19;14;12;23;24;22;1;18;21;15;11;3;10;5;9;20;17;6;7]; 
Lons = Lons(Order);
Lats = Lats(Order);


for iMarker=1:1:numel(Lats)
  m_text(Lons(iMarker),Lats(iMarker),Letters(iMarker),'fontsize',24,'fontweight','bold','color','k','horizontalalignment','center','verticalalignment','middle')
end; clear iMarker


%plot land
%%%%%%%%%%%%%%%%%%%%%%%
m_coast('color',[1,1,1].*0.3); %coasts at the bottom

%gridlines and done
%%%%%%%%%%%%%%%%%%%%%%%
m_grid('fontsize',14,'linestyle','--','color',[1,1,1].*.2);



% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % %% part B: find quarterly values for each box
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % 
% % % % % %what data do we want?
% % % % % Mode = 'h'; %cluster analysis
% % % % % Layer = 23; %lowest
% % % % % Stat =  3;  %median
% % % % % Range = [5,5]; %window size for average
% % % % % Variable = 'STT_A';
% % % % % 
% % % % % %load and merge quarterly files
% % % % % Store = struct();
% % % % % for iQ=1:1:4;
% % % % %   switch iQ; 
% % % % %     case 1; Q = 'son';  case 2; Q = 'djf'; case 3; Q = 'mam'; case 4; Q = 'jja';
% % % % %   end
% % % % %   LayerData = load(['out/',Mode,'_',Q,'_','b',num2str(Layer),'.mat']);
% % % % %   if iQ == 1; Store = LayerData.Results;
% % % % %   else        Store = cat_struct(Store,LayerData.Results,5);
% % % % %   end  
% % % % %   
% % % % %   if iQ == 1; LatScale = LayerData.Settings.Grid.Lat; LonScale = LayerData.Settings.Grid.Lon; end
% % % % % end; clear iQ LayerData Q
% % % % % 
% % % % % %pull out the maps we want
% % % % % Data = Store.(Variable);
% % % % % Data = squeeze(Data(:,:,1,Stat,:));
% % % % % 
% % % % % clear Mode Layer Stat Variable Store
% % % % % 
% % % % % Results = NaN(numel(Lats),4);
% % % % % for iMarker = 1:1:numel(Lats);
% % % % %   
% % % % %   InLatRange = inrange(LatScale,0.5.*Range(2)*[-1,1] + Lats(iMarker));
% % % % %   InLonRange = inrange(LonScale,0.5.*Range(1)*[-1,1] + Lons(iMarker));
% % % % %   
% % % % %   Box = Data(InLonRange,InLatRange,:);
% % % % % 
% % % % %   Results(iMarker,:) = nanmean(Box,[1,2]);
% % % % %   
% % % % % end; clear InLatRange InLonRange Box iMarker