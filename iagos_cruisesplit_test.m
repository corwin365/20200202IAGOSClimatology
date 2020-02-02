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
      Data = prep_iagos(Files{iFile},1./24./60./60.*Settings.dt);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% split data into "cruises", as suggested by Tim   
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
         
      %derivative of altitude
      dz = [diff(Data.baro_alt_AC),0];

      
      %produce a running windowed *sum* of this value
      dz = smooth(dz,Settings.Window).*Settings.Window;
      
      %take absolute value, as we're looking for discontinuities, not signs
      dz = abs(dz);
      
      %create a binary mask - 0 isgood, 1 is bad (during or near a height change)
      Track = zeros(size(dz)); Track(dz >= Settings.MaxDz) = 1;

      %and then divide the good bits into separate sequences seperated by
      %the bad bits. Remember that the plane will START AND END in a BAD
      %state, as it has to get to and from ground level.
      Starts = find(diff(Track) == -1); Starts = Starts(1:end-1);%last one is the end of the final descent
      Ends   = find(diff(Track) ==  1); 
      if numel(Starts) ~= numel(Ends); Ends = Ends(2:end); end %cross-check for sanity, as some have a bobble up and down at the start that triggers the end detector

      %these are the starts and ends of the cruises
      %split the data up
      csize = Ends-Starts+1;
      
      Cruises = NaN(numel(Starts),max(csize));
      for iCruise=1:1:numel(Starts);
        Cruises(iCruise,1:csize(iCruise)) = Starts(iCruise):1:Ends(iCruise);
      end
      clear iCruise Starts Ends Track dz csize
      

      
stop
    end
  end
end