clearvars


%find total duration of used cruises, to quote in the paper

TotalTime = 0;

for iYear=1994:1:2019;
  
  %load year
  ThisYear = [LocalDataDir,'/corwin/IAGOS_annual/merged_',num2str(iYear),'_sgolay900.mat'];
  load(ThisYear)
  
  %find duration of each cruise
  dt = range(Results.Time,2);
  TotalTime = TotalTime + nansum(dt(:));
  
end

%convert from days to hours
TotalTime = TotalTime.*24