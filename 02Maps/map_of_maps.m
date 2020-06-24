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


Names = {}; Lons = []; Lats = []; Oro = [];

Names{end+1} = 'Rockies';           Lons(end+1) = -110; Lats(end+1) = 40; Oro(end+1) = 1;
Names{end+1} = 'Iceland';           Lons(end+1) =  -18; Lats(end+1) = 65; Oro(end+1) = 1;
Names{end+1} = 'Newfoundland';      Lons(end+1) =  -60; Lats(end+1) = 48; Oro(end+1) = 0;
Names{end+1} = 'Greenland';         Lons(end+1) =  -47; Lats(end+1) = 64; Oro(end+1) = 1;
Names{end+1} = 'UK';                Lons(end+1) =   -2; Lats(end+1) = 54; Oro(end+1) = 1;
Names{end+1} = 'Iran';              Lons(end+1) =   49; Lats(end+1) = 35; Oro(end+1) = 1;
Names{end+1} = 'Altai';             Lons(end+1) =   90; Lats(end+1) = 52; Oro(end+1) = 1;
Names{end+1} = 'Sikhote Alin';      Lons(end+1) =  138; Lats(end+1) = 48; Oro(end+1) = 1;
Names{end+1} = 'Urals';             Lons(end+1) =   65; Lats(end+1) = 62; Oro(end+1) = 1;
Names{end+1} = 'Alps/Balkans';      Lons(end+1) =   15; Lats(end+1) = 46; Oro(end+1) = 1;
Names{end+1} = 'Canadian Plains';   Lons(end+1) = -110; Lats(end+1) = 61; Oro(end+1) = 0;
Names{end+1} = 'Siberia';           Lons(end+1) =  100; Lats(end+1) = 65; Oro(end+1) = 1;
Names{end+1} = 'Great Lakes';       Lons(end+1) =  -83; Lats(end+1) = 46.5; Oro(end+1) = 0;
Names{end+1} = 'North Atlantic';    Lons(end+1) =  -30; Lats(end+1) = 52; Oro(end+1) = 0;
Names{end+1} = 'Azores';            Lons(end+1) =  -25; Lats(end+1) = 37; Oro(end+1) = 0.5;
Names{end+1} = 'CMR Border';        Lons(end+1) =  120; Lats(end+1) = 55; Oro(end+1) = 1;
Names{end+1} = 'West Russia';       Lons(end+1) =   40; Lats(end+1) = 55; Oro(end+1) = 1;
Names{end+1} = 'Afghanistan';       Lons(end+1) =   69; Lats(end+1) = 38; Oro(end+1) = 1;


Letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

%order is computed in the small_multiples routine, and hardcode here for consistency
Order = 1:1:numel(Lons);%[1;14;13;16;4;18;2;17;8;6;7;10;12;15;3;9;5;11]
Lons = Lons(Order);
Lats = Lats(Order);


%place marker and draw 500km circle
for iMarker=1:1:numel(Lats)
  m_text(Lons(iMarker),Lats(iMarker),Letters(iMarker),'fontsize',24,'fontweight','bold','color','k','horizontalalignment','center','verticalalignment','middle')
  m_range_ring(Lons(iMarker),Lats(iMarker),500,'color','k','linewi',2)
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