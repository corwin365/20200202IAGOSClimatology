clearvars

%parameters to vary
Feed.Methods = {'h'};%{'g','h'};

Count = 1;

for iMethod =1:1:numel(Feed.Methods)
  for iBand=33%[1,5,7,31:1:35]
    for iDays=1:1:4;
      METHOD  = Feed.Methods{iMethod};
      PRSBAND = iBand;
      
      switch iDays
        case 1; DAYS = 'DJF';
        case 2; DAYS = 'MAM';
        case 3; DAYS = 'JJA';
        case 4; DAYS = 'SON';
      end
      
      OUTFILE = [METHOD,'_',DAYS,'_b',num2str(iBand),'_sgolay900'];
      
      
      JobName = OUTFILE;
      
      
      Text1 = ['#!/bin/bash\n## Name of the job\n#SBATCH --job-name=',JobName,'\n## Account to charge to\n#SBATCH --account=free\n\n'];
      Text2 = ['\n#SBATCH --time=240:00\n## Number of node required and tasks per node\n'];
      Text3 = ['#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=16\n\n#SBATCH --output=log.%%j.out\n#SBATCH --error=log.%%j.err\n#SBATCH --partition=batch-short'];
      
      
      %Load Matlab environment
      Text4 = ['\n\n\nmodule load matlab'];
      
      
      Text5 = ['\nmatlab -nodisplay -r "cd /home/f/cw785/Matlab/20200202IAGOSClimatology/02Maps;'];
      
      Text6 = ['METHOD = ''',METHOD,''';PRSBAND=',num2str(PRSBAND),';DAYS=''',DAYS,''';OUTFILE=''',OUTFILE,''';generate_map_data;exit"']
      
      fid = fopen(['job',sprintf('%04d',Count),'_wrapper.txt'],'wt');
      fprintf(fid, Text1);
      fprintf(fid, Text2);
      fprintf(fid, Text3);
      fprintf(fid, Text4);
      fprintf(fid, Text5);
      fprintf(fid, Text6);
      fclose(fid);
      
      Count = Count+1;
    end
  end
end
  
%generate file to fire the whole lot off
fid = fopen(['fire_wrappers.sh'],'wt');
%  for i=1:1:Count-1;fprintf(fid,['sbatch --begin=now+7hour job',sprintf('%04d',i),'.txt\n']);end
for i=1:1:Count-1;
  fprintf(fid,['sbatch  job',sprintf('%04d',i),'_wrapper.txt\n']);
  fprintf(fid,['rm job',sprintf('%04d',i),'_wrapper.txt\n']);
end
fprintf(fid,['rm fire_wrappers.sh\n']);
fclose(fid);

disp(['Written ',num2str(Count),' files (probably)'])

