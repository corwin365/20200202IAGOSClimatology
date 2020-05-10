clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%where is the data?
Settings.DataDir = [LocalDataDir,'/IAGOS/Timeseries']

%what time step to interpolate to?
Settings.dt = 5; %second

%how is a 'cruise' defined?
  %first variable is height change required
  %second variable is how long this change is over, in the time units specified above
Settings.MaxDz  = 100;%m
Settings.Window = 15*60./Settings.dt;  %fifteen minutes

%detrending window size
Settings.Detrend = 30*60./Settings.dt;

%smoothing size
Settings.Smooth = 5*60./Settings.dt;



%mapping of final results
Settings.Lon = -130:.5:30;
Settings.Lat = 10:.5:90;


Settings.Years = [1994,2006];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Settings.Window = round(Settings.Window);
if ~isodd(Settings.Window); Settings.Window = Settings.Window+1; end

Settings.Detrend = round(Settings.Detrend);
if ~isodd(Settings.Detrend); Settings.Detrend = Settings.Detrend+1; end

Settings.Smooth = round(Settings.Smooth);
if ~isodd(Settings.Smooth); Settings.Smooth = Settings.Smooth+1; end


[xi,yi] = meshgrid(Settings.Lon,Settings.Lat);

Results.Sum = zeros([size(xi),numel(Settings.Years)]);
Results.Count = Results.Sum;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iYear=1:1:numel(Settings.Years)

  for iMonth=1:2;%:1:12;
    
    Files = wildcardsearch([Settings.DataDir,...
                           '/',sprintf('%02d',Settings.Years(iYear)), ...
                           '/',sprintf('%02d',iMonth),'/'],'*.nc');

    textprogressbar(['Processing ',sprintf('%04d',Settings.Years(iYear)),'/',sprintf('%02d',iMonth),' '])
    for iFile=1:1:numel(Files);
      textprogressbar(iFile./numel(Files).*100)
      

      try
      %load file, including unit conversions
      Data = prep_iagos(Files{iFile}, ...
                        'SamplingRate',1./24./60./60.*Settings.dt, ...
                        'CruiseDz',Settings.MaxDz, 'CruiseWindow',Settings.Window);

      %loop over cruises
      for iCruise = 1:1:size(Data.Cruises,1)
       
        %identify cruise
        Cruise = Data.Cruises(iCruise,:);
        Cruise = Cruise(~isnan(Cruise));
       
        %extract data and 
        Var = Data.air_temp_AC(Cruise);
        Lon = Data.lon(Cruise); 
        Lat = Data.lat(Cruise); 
        
        %smooth out small scales
        Var = smooth(Var,Settings.Smooth);
        
        
        %detrend large scales
        BG = smooth(Var,Settings.Detrend);
        Var = Var - BG;
        
        %s-transform
        ST = nph_ndst(Var);
        
        %put amplitude onto map
        ZG  = bin2matN(2,Lon,Lat,ST.A,xi,yi);
        ZGc = bin2matN(2,Lon,Lat,ones(size(ST.A)),xi,yi,'@sum'); 
        
        ZG(isnan(ZG)) = 0;
        ZGc(isnan(ZGc)) = 0;
        
        
        Results.Sum(  :,:,iYear) = Results.Sum(  :,:,iYear) + ZG;
        Results.Count(:,:,iYear) = Results.Count(:,:,iYear) + ZGc;
        

       
      end
      
      catch; end
    end
    textprogressbar('!')
   
    save('out2.mat','Results','xi','yi','Settings')
  end
end