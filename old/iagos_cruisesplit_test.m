clearvars

%% settings

%where is the data?
Settings.DataDir = [LocalDataDir,'IAGOS/TimeSeries'];

%what time step to interpolate to?
Settings.dt = 5; %second

%how is a 'cruise' defined?
  %first variable is height change required
  %second variable is how long this change is over, in the time units specified above
Settings.MaxDz  = 100;%m
Settings.Window = 15*60./Settings.dt;  %fifteen minutes


%% preprocessing
Settings.Window = round(Settings.Window);
if ~isodd(Settings.Window); Settings.Window = Settings.Window+1; end


%% processing


for iYear=2006

  for iMonth=1:1:12;
    
    Files = wildcardsearch([Settings.DataDir,...
                           '/',sprintf('%02d',iYear), ...
                           '/',sprintf('%02d',iMonth),'/'],'*.nc');

    for iFile=1:1:numel(Files);

      
      %load file, including unit conversions
      Data = prep_iagos(Files{iFile}, ...
                        'SamplingRate',1./24./60./60.*Settings.dt, ...
                        'CruiseDz',Settings.MaxDz, 'CruiseWindow',Settings.Window);

   

    end
  end
end