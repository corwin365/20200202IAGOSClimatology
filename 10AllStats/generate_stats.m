clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate statistics of GW properties measured by IAGOS
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/Jun/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Settings.DataDir    = [LocalDataDir ,'/corwin/IAGOS_annual/'];
Settings.TimeScale  = datenum(1994,8,1):1:datenum(2019,12,31);
% Settings.TimeScale  = datenum(2006,1,1):1:datenum(2006,12,31);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load all data into one pile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Store.OldYear = 0;

for Time = min(Settings.TimeScale):1:max(Settings.TimeScale);

  [yy,~,~] = datevec(Time);
  if Store.OldYear ~= yy

    %load data
    FileName = [Settings.DataDir,'/merged_',num2str(yy),'_sgolay900.mat'];
    if ~exist(FileName,'file'); continue; end
    File = load(FileName); File = File.Results;
    
    %create an array for unique flight IDs
    File.NU = repmat(1:1:size(File.Lat,1),size(File.Lat,2),1)';
    
    %require non-nan data
    Good = find(~isnan(File.Lat + File.Lon + File.STT_A));

%     %require wavelength greater than 80km
%     Good = intersect(Good,find(1./File.STT_k > 25));
    
    %glue to our arrays
    if ~isfield(Store,'Lat'); Store = reduce_struct(File,Good);
    else;                     Store = cat_struct(Store,reduce_struct(File,Good),1,{'OldYear'});
    end
    Store.OldYear = yy;
    clear InLatRange InLonRange InRange File InPrsRange Good
    disp(['Loaded ',num2str(yy)]);
  end
              
end; clear TimeRange yy FileName Time
disp('Data loaded')

Store.Lambda = (1./Store.STT_k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prepare figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf
set(gcf,'color','w')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% temperature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,1)

x = linspace(190,270,100);
h = histogram(Store.T,x);

xlim([min(x),max(x)])
ylabel('Count'); xlabel('Temperature [K]')
title('(a) Temperature');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% zonal wind speed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,2)

x = linspace(-50,100,100);
h = histogram(Store.U,x);

xlim([min(x),max(x)])
ylabel('Count'); xlabel('Wind Speed [m/s]')
title('(b) Zonal Wind Speed');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% meridional wind speed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,3)

x = linspace(-75,75,100);
h = histogram(Store.V,x);

xlim([min(x),max(x)])
ylabel('Count'); xlabel('Wind Speed [m/s]')
title('(c) Meridional Wind Speed');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% amplitude histogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,5)

x = logspace(-1.5,1,100);
h = histogram(Store.STT_A,x);

xlim([min(x),max(x)])
set(gca,'xscale','log')
ylabel('Count'); xlabel('Amplitude [K]')
title('(e) Amplitudes');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% wavenumber histogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,6)

x = logspace(log10(10),log10(1000),100);
h = histogram(Store.Lambda(Store.Lat < 50),x);
% plot(x,y)

xlim([min(x),max(x)])
set(gca,'xscale','log')
ylabel('Count'); xlabel('Wavelength [km]')
title('(f) Wavelengths');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pressure histogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,4)

x = linspace(150,350,50);
h = histogram(Store.Prs,x);

set(gca,'yscale','log');% ylim([10^2 10^6.5])
ylabel('Count'); xlabel('Pressure [hPa]')
title('(d) Pressure Levels');
