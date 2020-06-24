clearvars -except YEAR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%general setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%where is the data?
Settings.DataDir = [LocalDataDir,'/IAGOS/Timeseries'];
Settings.OutDir  = [LocalDataDir,'/corwin/IAGOS_st/'];

%dates to loop over. A separate file will be produced for each day.
Settings.TimeScale = datenum(YEAR,1,1):1:datenum(YEAR,12,31);

%cruise identification
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%what time step to interpolate to?
Settings.dt = 5; %second

%how is a 'cruise' defined?
%first variable is height change required
%second variable is how long this change is over, in the time units specified above
Settings.MaxDz  = 100;%m
Settings.Window = 15*60./Settings.dt;  %fifteen minutes

%minimum length (km) of cruise
Settings.MinCruiseLength = 1000;

%maximum size of any discontinuities before complete discard of dataset
Settings.MaxDiscontinuity = 200;% km

%spectral analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%distance spacing
Settings.SA.dx = 1; %km

%maximum gap size
Settings.SA.MaxSpaceGap = 20./Settings.SA.dx; %km

%low-pass filter size
Settings.SA.Detrend = 500./Settings.SA.dx;

%minimum wavelengths
Settings.MinLambda = 20; %km

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% variables to retain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%geolocation
Settings.Vars.Geo.In  = {'lat','lon','Time'};
Settings.Vars.Geo.Out = {'Lat','Lon','Time'};

%s-transform - temperature
Settings.Vars.STT.In  = {'IN','A','F1','EdgeMask'};
Settings.Vars.STT.Out = {'Tprime','STT_A','STT_k','STT_EdgeMask'};

%s-transform - U
Settings.Vars.STU.In  = {}%{'IN','A','F1','EdgeMask'};
Settings.Vars.STU.Out = {}%{'Uprime','STU_A','STU_k','STU_EdgeMask'};

%s-transform - V
Settings.Vars.STV.In  = {}%{'IN','A','F1','EdgeMask'};
Settings.Vars.STV.Out = {}%{'Vprime','STV_A','STV_k','STV_EdgeMask'};

%metadata
Settings.Vars.Meta.In  = {'air_press_AC','baro_alt_AC','zon_wind_AC','mer_wind_AC','air_temp_AC'};
Settings.Vars.Meta.Out = {'Prs','Z','U','V','T'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Settings.Window = round(Settings.Window);
if ~isodd(Settings.Window); Settings.Window = Settings.Window+1; end

Settings.SA.Detrend = round(Settings.SA.Detrend);
if ~isodd(Settings.SA.Detrend); Settings.SA.Detrend = Settings.SA.Detrend+1; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iDay =1:1:numel(Settings.TimeScale)
  
  OutFile = [Settings.OutDir,'/IAGOS_ST_',num2str(Settings.TimeScale(iDay)),'_sgolay.mat'];
  
  if exist(OutFile); 
    %check when file was last modified 
    file = dir(OutFile);
    %last complete rerun was on the date below, so ignore any date after that...
    if datenum(file.date) > datenum(2020,5,27,23,30,00);     
      disp([datestr(Settings.TimeScale(iDay)),' already done'])
      continue; 
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %create results arrays for day
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  X = NaN(1,1);%will expand out as needed
  Results = struct();
  Vars = [Settings.Vars.Geo.Out, ...           
          Settings.Vars.Meta.Out, ...
          Settings.Vars.STT.Out, ...
          Settings.Vars.STU.Out, ...
          Settings.Vars.STV.Out,];
  for iVar=1:1:numel(Vars);
    Results.(Vars{iVar}) = X;
  end
  clear iVar Vars X
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %load day's data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  [yy,mm,dd] = datevec(Settings.TimeScale(iDay));
  disp(datestr(Settings.TimeScale(iDay)))
  
  
  
  Files = wildcardsearch([Settings.DataDir,'/iagos_',sprintf('%04d',yy)],...
                          ['*',sprintf('%04d',yy),sprintf('%02d',mm),sprintf('%02d',dd),'*.nc']);
  for iFile=1:1:numel(Files);
    
    
try
      %load file, including unit conversions
      
      %in this step we interpolate to time to identify the cruises. space
      %intepolation will be done instead in the loop below to get accurate
      %spatial wavelengths
      try
        Data = prep_iagos(Files{iFile}, ...
                          'SamplingRate',1./24./60./60.*Settings.dt, ...
                          'CruiseDz',Settings.MaxDz, 'CruiseWindow',Settings.Window, ...
                          'ApplyFlags',true);
      catch;
        disp(['Error in file ',Files{iFile}])
        continue
      end

      %loop over cruises
      for iCruise = 1:1:size(Data.Cruises,1)
 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %find cruise info and regularly grid
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        %identify cruise
        Cruise = Data.Cruises(iCruise,:);
        Cruise = Cruise(~isnan(Cruise));
                
        %extract geolocation
        Lat = Data.lat(Cruise);    
        Lon = Data.lon(Cruise);     
        
        %get distances for regular interpolation
        x = [Lat;Lon]'; y = circshift(x,1,1);
        dx = nph_haversine(x,y); dx(1) = dx(2);
        dxS = dx; for iX=2:1:numel(dxS); dxS(iX) = dxS(iX-1)+dxS(iX); end  
        dx2 = 0:Settings.SA.dx:max(dxS);
        
        %check all distance elements are unique (drop any duplicates below)
        if find(abs(dx) <= 1e-3); [~,Unique] = unique(dxS,'stable'); %within a metre is as good as zero on these scales, and fixes some edge cases where 0 doesn't work
        else;                     Unique = 1:1:numel(dxS);
        end

        %check all distance points are non-NaN
        Good = find(~isnan(dxS(Unique)));
        Unique = Unique(Good);
        clear Good

        %check we still have enough data to be useful
        if numel(Unique) < 10; continue; end        
        
        %check the cruise is long enough
        if max(dxS) < Settings.MinCruiseLength; continue; end

        %interpolate all variables to fixed grid
        Vars = fieldnames(Data);
        Regular = struct();
        for iVar=1:1:numel(Vars)
          if strcmp(Vars{iVar},'MetaData');     continue; end
          if strcmp(Vars{iVar},'OriginalTime'); continue; end
          if strcmp(Vars{iVar},'Cruises');      continue; end
          
          Var = Data.(Vars{iVar});
          if sum(~isnan(Var(Cruise(Unique)))) <2 ; continue; end
          
          Regular.(Vars{iVar}) = interp1gap(dxS(Unique),Var(Cruise(Unique)),dx2,Settings.SA.MaxSpaceGap);
          
          %scale pressures
          if strcmp(Vars{iVar},'Prs'); Regular.(Vars{iVar}) = Regular.(Vars{iVar})./100; end
          
        end
        Regular.dxS = dx2;
        clear iVar Vars dx dxS dx2 x y
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %detrend and s-transform T, U and V
        %
        %do each variable separately, so we can check if they exist
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if ~isfield(Regular,'air_temp_AC'); continue; end        
        
        for iVar=1%1:1:3;
        
          try
            %get the variable
            switch iVar
              case 1; Var = Regular.air_temp_AC;
              case 2; Var = Regular.zon_wind_AC;
              case 3; Var = Regular.mer_wind_AC;
            end
          catch
            %variable must not exist
            continue
          end
        
          %trim any NaNs at the end - these can cause artificially large wave
          %modes to be fitted due to the nanfilling step later
          %The NaNs will be put back at the end so that the data correspond
          %to the other series in length
          NaNTrim.NElements = numel(Var);
          NaNTrim.Keep = find(~isnan(Var));
          Var = Var(NaNTrim.Keep);
          
          
          %detrend large scales
% %           Var = Var-smoothdata(Var,'gaussian',Settings.SA.Detrend);
          FrameSize = Settings.SA.Detrend; if ~isodd(FrameSize); FrameSize = FrameSize+1; end
          Var = Var-sgolayfilt(double(Var),2,Settings.SA.Detrend);
          
          %remove if there are large discontinuities
          Disco = diff(find(~isnan(Var))).*Settings.SA.dx;
          if max(Disco) > Settings.MaxDiscontinuity; continue; end
          
          %fill all the gaps for the s-transform (undone after)
          Bad = find(isnan(Var));
          Var = inpaint_nans(double(Var));
          
          %s-transform
          %scales: only use those in the valid range of wavelengths after filtering
          LambdaRange = [Settings.MinLambda,500];%Settings.SA.Detrend.*Settings.SA.dx];
          Len = numel(Var).*Settings.SA.dx;
          Lambdas = min(LambdaRange):Settings.SA.dx:max(LambdaRange);
          Scales = unique(Len./Lambdas);
          
          %do s-transform
          ST = nph_ndst(Var,Scales,Settings.SA.dx); %simple 1dst
                    
          %put the gaps back
          ST.In(Bad) = NaN;
          ST.A( Bad) = NaN;
          ST.F1(Bad) = NaN;
          
          %mask anything edge-truncated
          CutOff = abs((1:1:numel(ST.F1))-numel(ST.F1)./2); CutOff = CutOff - CutOff(1);
          Over   = 1./ST.F1 + 2.*CutOff;
          ST.EdgeMask = zeros(size(ST.In));
          ST.EdgeMask(Over > 0) = 1;
                   
          %finally, put back the end-NaNs
          A2 = NaN(NaNTrim.NElements,1); A2(NaNTrim.Keep) = ST.A;  ST.A  = A2;
          k2 = NaN(NaNTrim.NElements,1); k2(NaNTrim.Keep) = ST.F1; ST.F1 = k2;
          i2 = NaN(NaNTrim.NElements,1); i2(NaNTrim.Keep) = ST.IN; ST.IN = i2; 
          e2 = NaN(NaNTrim.NElements,1); e2(NaNTrim.Keep) = ST.EdgeMask; ST.EdgeMask = e2;

          %and store
          switch iVar
            case 1; STT = ST;
            case 2; STU = ST;
            case 3; STV = ST;
          end
          
          %tidy up
          clear LambdaRange Len Scales A2 k2 ST Lambdas Bad Var Disco NaNTrim Over

        end; clear iVar
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %store the outputs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for VarType = {'Meta','Geo','STT','STU','STV'};
          
          if any(strcmp(VarType{1},{'STT','STU','STV'})) == 1 & ~exist(VarType{1},'var') 
            %variable does not exist. This is often the case for valid 
            %reasons - e.g. a good T time series is available but the U
            %or V data were discarded
            continue;
          end
          
          for iVar=1:1:numel(Settings.Vars.(VarType{1}).Out)
           
            %pull results array out of struct
            R = Results.(Settings.Vars.(VarType{1}).Out{iVar});
            
            %get data to store
            switch VarType{1};
              case 'STT'; O = STT;
              case 'STU'; O = STU;
              case 'STV'; O = STV;
              otherwise;  O = Regular;
            end
              
            %if the field doesn't exist, set it all to NaNs
            if isfield(O,Settings.Vars.(VarType{1}).In{iVar});
              O = O.(Settings.Vars.(VarType{1}).In{iVar});
            else
              O = NaN(size(Regular.UTC_time)); %should always be there
            end
            
            %store data, and fill any end-gaps with NaNs
            oldsize = size(R,2);
            R(iCruise,1:numel(O)) = O;
            if oldsize > numel(O);
              R(iCruise,numel(O)+1:oldsize) = NaN;
            end
            
            %put results back into struct
            Results.(Settings.Vars.(VarType{1}).Out{iVar}) = R;
          end; clear R O iVar oldsize
          
        end
      end; clear iCruise
      
      

 catch; 
   disp(['Error on ',datestr(Settings.TimeScale(iDay))])
 end
  end; clear iFile
  
% % % % % %   %because STU, STT and STV aren't generated every loop, it's possible
% % % % % %   %that they may have ended up as shorter time series
% % % % % %   %check this, and NaN-pad them if so to help logic of later routines
% % % % % %   Longest = max([size(Results.STT_A,2),size(Results.STU_A,2),size(Results.STV_A,2)]);
% % % % % %   
% % % % % %   
% % % % % %   if isfield(Results,'STU_A')
% % % % % %     if size(Results.STU_A,2) < Longest;
% % % % % %       Extra = NaN(size(Results.STU_A,1), ...
% % % % % %                   Longest - size(Results.STU_A,2));
% % % % % %       Results.STU_A        = cat(2,Results.STU_A, Extra);
% % % % % %       Results.STU_k        = cat(2,Results.STU_k, Extra);
% % % % % %       Results.STU_EdgeMask = cat(2,Results.STU_EdgeMask, Extra);      
% % % % % %       Results.Uprime       = cat(2,Results.Uprime,Extra);
% % % % % %     end
% % % % % %   end
% % % % % %   if isfield(Results,'STV_A')
% % % % % %     if size(Results.STV_A,2) < Longest;
% % % % % %       Extra = NaN(size(Results.STV_A,1), ...
% % % % % %                   Longest - size(Results.STV_A,2));
% % % % % %       Results.STV_A        = cat(2,Results.STV_A, Extra);
% % % % % %       Results.STV_k        = cat(2,Results.STV_k, Extra);
% % % % % %       Results.STV_EdgeMask = cat(2,Results.STV_EdgeMask, Extra);  
% % % % % %       Results.Vprime       = cat(2,Results.Vprime,Extra);
% % % % % %     end
% % % % % %   end
% % % % % %   if isfield(Results,'STT_A')
% % % % % %     if size(Results.STT_A,2) < Longest;
% % % % % %       Extra = NaN(size(Results.STT_A,1), ...
% % % % % %                   Longest - size(Results.STT_A,2));
% % % % % %       Results.STT_A        = cat(2,Results.STT_A, Extra);
% % % % % %       Results.STT_k        = cat(2,Results.STT_k, Extra);
% % % % % %       Results.STT_EdgeMask = cat(2,Results.STT_EdgeMask, Extra);  
% % % % % %       Results.Tprime       = cat(2,Results.Tprime,Extra);
% % % % % %     end
% % % % % %   end  
% % % % % %   clear Extra Longest
% % % % % %   

  %finally, store the data for the day
  if nansum(Results.Lon(:)) ~= 0;  save(OutFile,'Results','Settings'); end

end; clear iDay
