clearvars


%find total duration of used cruises, to quote in the paper

TotalTime = 0;

for iYear=1994:1:2020;
  
  ThisYear = [LocalDataDir,'/corwin/IAGOS_annual/merged_',num2str(iYear),'_sgolay900.mat'];
  load(ThisYear)
  
  stop
end