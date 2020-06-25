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
Generate = 1; %only generate new data if needed, it is *very* slow

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

%order is hardcoded to be the same as the small multiples annualised plots
Order = [1;14;3;5;13;6;10;18;4;11;16;2;15;7;8;17;12;9];
Lons  = Lons(Order);
Lats  = Lats(Order);
Names = Names(Order);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate the time series (if needed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Generate == 1;
  
  for iSeries = 1%:1:numel(Names);
    disp(['Processing box over ',Names{iSeries}])
    
    %generate the time series
    func_generate_composite_time_series_lt(urlencode(Names{iSeries}),Lons(iSeries),Lats(iSeries),Range, PrsRange,TimeScale, 'STT_A', 'nanmean', TimeWindow)

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


