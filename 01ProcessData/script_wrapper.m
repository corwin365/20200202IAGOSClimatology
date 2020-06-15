clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% balena operation settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%name of the programme to be used
ToCall.ProgName = 'combine_data_into_annual_files';

%directory it will run in
ToCall.WorkingDirectory = '/home/f/cw785/Matlab/20200202IAGOSClimatology/01ProcessData';

Queue = 'batch-devel'  ;%options are batch-short, batch-sky, batch-all, batch-devel
RunTime = 360; %this will be overriden with 15 min if 'batch-devel' queue is used, and only one job fired


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate individual function calls
%
%this part will be different for every job. Aim is to create a
%cell array with each entry representing an individual job
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Calls = {};
JobNames = {};

Count = 0;
for YEAR = 1994:1:2020
    Count = Count+1;

    ThisCall = ['YEAR=',num2str(YEAR),';', ...
               'iagos_cruises_st_fixedDX;add_tropopause;'];

    JobNames{Count} = ['IAGOS_',sprintf('%04d',YEAR)];
    Calls{Count} = ThisCall;
end
clearvars -except Calls ToCall Queue RunTime JobNames





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate individual job scripts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%special settings for certain queues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%number of cores
switch Queue
  case 'batch-sky';   Cores = 24;
  otherwise;    Cores = 16;
end

%quality of service
switch Queue
  case 'batch-devel'; QOS = '#SBATCH --qos=devel\n';
  otherwise;    QOS = ''; 
end

%max runtime
switch Queue
  case 'batch-devel'; RunTime = 15;
  otherwise; %already set, don't override
end

%max jobs
switch Queue
  case 'batch-devel'; MaxJobs = 1;
  otherwise;          MaxJobs = numel(Calls);
end

%ok, make the scripts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Count = 1;
for iScript=1:1:MaxJobs
  
  %general setup for Matlab balena jobs the random pause at the end 
  %is to spread out parpool calls if used in the scripts, as multiple calls
  %in rapid succession can cause scripts to crash when creating the pool 
 Start = ['#!/bin/bash\n', ...
         '#SBATCH --job-name=',JobNames{iScript},'\n', ...
         '#SBATCH --account=free\n', ...
         '#SBATCH --time=',num2str(RunTime),':00\n', ...
         '#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=',num2str(Cores),'\n', ...
         QOS, ...
         '#SBATCH --output=log.%%j.out\n', ...
         '#SBATCH --error=log.%%j.err\n', ...
         '#SBATCH --partition=',Queue,'\n', ...
         'module load matlab\n',...
         'matlab -nodisplay -r "cd ',ToCall.WorkingDirectory,';', ...
         'pause(',num2str(randi(30)),');'];
       End = ';exit;"';

 
  Commands = [Calls{iScript},';',ToCall.ProgName];
  
  %create scripts
  fid = fopen(['job',sprintf('%04d',Count),'_wrapper.txt'],'wt');
  fprintf(fid, Start);fprintf(fid, Commands);fprintf(fid, End);
  fclose(fid);
  Count = Count+1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate master script to call the above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%generate file to fire the whole lot off
fid = fopen(['fire_wrappers.sh'],'wt');
for i=1:1:Count-1;
  fprintf(fid,['sbatch job',sprintf('%04d',i),'_wrapper.txt\n']);
  fprintf(fid,['rm job',sprintf('%04d',i),'_wrapper.txt\n']);
end
fprintf(fid,['rm fire_wrappers.sh\n']);
fclose(fid);

disp(['Written ',num2str(Count),' files (probably)'])

