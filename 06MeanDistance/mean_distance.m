clearvars


%find basic statistics on original data spacing, to quote in the paper

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%general setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%where is the data?
Settings.DataDir = [LocalDataDir,'/IAGOS/Timeseries'];

%dates to loop over. A separate file will be produced for each day.
Settings.TimeScale = datenum(1994,1,1):1:datenum(2019,12,31);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AllStore =  NaN(2e9,1); %arbitrarily big!
Count = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iDay =1:1:numel(Settings.TimeScale)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %load day's data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  [yy,mm,dd] = datevec(Settings.TimeScale(iDay));
  disp(datestr(Settings.TimeScale(iDay)))
  
   
  Files = wildcardsearch([Settings.DataDir,'/iagos_',sprintf('%04d',yy)],...
                          ['*',sprintf('%04d',yy),sprintf('%02d',mm),sprintf('%02d',dd),'*.nc']);
  for iFile=1:1:numel(Files);

      %load file, including unit conversions
      
      %in this step we interpolate to time to identify the cruises. space
      %intepolation will be done instead in the loop below to get accurate
      %spatial wavelengths
      Data = rCDF(Files{iFile});
      
      %get distances for regular interpolation
      x = [Data.lat,Data.lon]; y = circshift(x,1,1);
      dx = nph_haversine(x,y); dx(1) = dx(2);

      %store
      idx = Count+1:Count+numel(dx);
      AllStore(idx) = dx;
      Count = Count+numel(idx);
      
      
      clear x dx idx y
      
      clear Data
  end; clear iFile Files
  
  clear yy mm dd 
end; clear iDay
 

Spacing = nanmean(AllStore)
StDev   = nanstd(AllStore)