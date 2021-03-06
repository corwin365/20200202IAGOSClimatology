clearvars


disp('does not work')
stop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%find tropopause in ERA5 data
%%
%%
%%approximates the method of Reichler (2003, GRL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%file and time handling
Settings.DataDir = [LocalDataDir,'/ERA5'];
Settings.OutFile = 'era5_tropopause.mat';
Settings.TimeScale = datenum(2010,1,1):datenum(2010,1,5);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C.R  = 287;
C.cp = 1;
C.g  = 9.81;

C.K = C.R ./ C.cp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% processing (output files created inside loop on first iteration)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iDay=1:1:numel(Settings.TimeScale)
  
  
  %find and load the ERA5 data for this day
  [yy,~,~] = datevec(Settings.TimeScale(iDay));
  dd = date2doy(Settings.TimeScale(iDay));
  File = [Settings.DataDir,'/',sprintf('%04d',yy),'/',...
          'era5_',sprintf('%04d',yy),'d',sprintf('%03d',dd),'.nc'];
  if ~exist(File); continue; end
  Data = rCDF(File);
  Prs = ecmwf_prs_v2([],137);
  clear yy dd File
  
  %store geoloc for later
  Lat = Data.latitude;
  Lon = Data.longitude;
  
  %extract temperature
  T = Data.t;
  clear Data
  
  %make the data go in the right direction
  [~,idx] = sort(Prs,'desc');
  Prs = Prs(idx);
  T = T(:,:,idx,:);
  clear idx
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% algorithm step 1: lapse rate 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %find half-level pressures
  Ph = 0.5 .*(Prs + circshift(Prs,1,1));
  Ph = Ph(2:end); %can't use due to circshift
  
  %find lapse rate at each half-level
  P = permute(repmat(permute(Prs,[2,3,4,1]),size(T,1),size(T,2),size(T,4),1),[1,2,4,3]);
  Term1 = diff(T,1,3)./diff(P,1,3);
  
  Term2 = (P + circshift(P,1,3)) ./ (T + circshift(T,1,3));
  Term2 = Term2(:,:,2:end,:);
  
  Term3 = C.K .* C.g ./ C.R;
  
  Gamma = Term1 .* Term2 .* Term3;
  stop
  clear Term1 Term2 Term3 T Prs P

  
  stop
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% algorithm step 2: find tropopause
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  %create an array to store our tropopause levels
  sz = size(Gamma);
  Tropopause = NaN(sz([1,2,4]));
  clear sz
  
  
  %loop over half-levels
  textprogressbar('Finding tropopause ')
  for iLevel = 1:1:136
    textprogressbar(iLevel./136.*100)
    %if we've already found all the tropopauses, continue
    if sum(isnan(Tropopause(:))) == 0; continue; end
    %if half-pressure > 550hPa, skip
    if Ph(iLevel) > 550; continue; end
    %if half-pressure < 75, skip
    if Ph(iLevel) < 75; continue; end

    %check if Gamma is less than -2
    idx = find(Gamma(:,:,iLevel,:) < -2);
    
    %remove any where we already found the tropopause
    if nansum(Tropopause) ~=0; stop; end
    
    if numel(idx) == 0; continue; end %none at this level
    
    %for each element where the above criterion is met, check if the layer
    %2km higher also meets it
    
    %find which level is 2km above
    Z = p2h(Ph(iLevel));
    jLevel = closest(p2h(Ph),Z+2);
    Range = sort([iLevel,jLevel],'ascend'); 
    Range = Range(1):1:Range(2);
    clear Z jLevel
    
    %pull out this range for each of the elements of interest
    G2 = permute(Gamma(:,:,Range,:),[1,2,4,3]);
    sz = size(G2);
    G2 = reshape(G2,sz(1)*sz(2)*sz(3),sz(4));
    G2 = G2(idx,:);
    clear Range sz
    
    %find all the columns where the criterion remains met for 2km above
    StillMet = min(G2,[],2);
    Good = find(StillMet < -2);
    idx = idx(Good);
% %     stop
% %     
% %     %ok. all of these meet the criterion
% %     %now, for each one find the actual tropopause via interpolation
% %     stop    

    
    %yay :-) store, and remove these columns from the lapse rate data
    Tropopause(idx) = Ph(iLevel);
    
    
    
  end; clear iLevel
  textprogressbar(1); textprogressbar('!')
  stop
    
  
end


