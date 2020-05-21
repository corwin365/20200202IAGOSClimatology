clearvars


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% merge IAGOS-ST daily files into annual files, to reduce overhead in data loading
%
%Corwin Wright, c.wright@bath.ac.uk
%2020/05/20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AnnSettings.InDir  = [LocalDataDir,'/corwin/IAGOS_st/'];
AnnSettings.OutDir = [LocalDataDir,'/corwin/IAGOS_annual/'];
AnnSettings.Years  = 2004:2011;%2016:-1:2014; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iYear=1:1:numel(AnnSettings.Years)
  
  disp(['Merging ',num2str(AnnSettings.Years(iYear))])
  

  %loop over days in the year, loading files
  %do this in steps of 40 days, to keep the array sizes mangeable, then
  %merge at the end
  textprogressbar('--> Processing days ')
  First = 1;
  DataThisYear = 0;
  for iDay = 1:1:365 %boohoo to leap days
    textprogressbar(iDay./365.*100)
    
    if mod(iDay,40) == 1; 
      %tidy up previous loop
      if First ~= 1 && DataThisYear ~=0;
        Store.(['r',num2str(iDay)]) = Results;
      end      
      clear Results Count
      First = 0;
    end
    
    %find the day of data
    DayNumber = datenum(AnnSettings.Years(iYear),1,iDay);
    
    %load the file
    FileName = [AnnSettings.InDir,'/IAGOS_ST_',num2str(DayNumber),'_v4.mat'];
    clear DayNumber
    if ~exist(FileName,'file'); clear FileName; continue; end
    
    DataThisYear = 1;
    Data = load(FileName);
    
    %get list of variables
    Vars = fieldnames(Data.Results);
    
    if ~exist('Results','var')
      %if this is the first day, create a results array that is large,
      
      for iVar=1:1:numel(Vars)
        Results.(Vars{iVar}) = NaN(400,5000); %400 records of max length 5000
        Count.(  Vars{iVar}) = 0; %number of records so far
      end
      
      %and also extract the Settings struct, so we can re-store it
      Settings = Data.Settings;
        
    end; clear iVar
    
    %check how big the data is, and store
    szNew = size(Data.Results.Lon);
    
    for iVar=1:1:numel(Vars)
      V = Results.(Vars{iVar});
      C = Count.(  Vars{iVar});
      
      V(C+1:C+szNew(1), 1:szNew(2)) = Data.Results.(Vars{iVar});
      Results.(Vars{iVar}) = V;  
      Count.(  Vars{iVar}) = C+szNew(1);
    end
    
    clear iVar szNew C V Vars
  end; 
  if exist('Results'); Store.(['r',num2str(iDay)]) = Results; end
  clear iDay
  clear Results Count 
  textprogressbar('!')


  %merge the 50-day chunks
  Chunks = fieldnames(Store);
  Results = Store.(Chunks{1});
  for iChunk=2:1:numel(Chunks)
    Vars = fieldnames(Store.(Chunks{1}));
    for iVar=1:1:numel(Vars)
      Results.(Vars{iVar}) = cat(1,Results.(Vars{iVar}),Store.(Chunks{iChunk}).(Vars{iVar}));
      
    end
  end
  clear Chunks iChunk Vars iVar Store

  
  %drop empty rows
  GoodX = find(nansum(Results.Lon,1) ~=0);
  GoodY = find(nansum(Results.Lon,2) ~=0);
  Vars = fieldnames(Results);
  for iVar=1:1:numel(Vars)
    V = Results.(Vars{iVar});
    V = V(GoodY,:);
    V = V(:,GoodX);
    Results.(Vars{iVar}) = V;
  end
  clear GoodX GoodY Vars iVar V

  
  %save the big annual pile of data
  save([AnnSettings.OutDir,'/merged_',num2str(AnnSettings.Years(iYear)),'.mat'], ...
       'Results','Settings');

     
     
end; clear iYear