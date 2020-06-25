clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% map indicating locations discussed in text
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


Letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

%order is computed in the small_multiples routine, and hardcoded here by hand
Order = [1;14;3;5;13;6;10;18;4;11;16;2;15;7;8;17;12;9]
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



