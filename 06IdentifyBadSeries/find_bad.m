clearvars

%find profiles that are bad, then flag them for fixing
%the bug or oversight that let each of these through the processing chain
%has been fixed, and the relevant dates re-run. Therefore, this file is
%hisootircal and not used in the main analysis chain for the paper.

Dates = [];
for iYear=1994:1:2019
  disp(iYear)
  
  %load the year
  load([LocalDataDir,'/corwin/IAGOS_annual/merged_',num2str(iYear),'.mat'])
  
% %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   %%
% %   %problem 1: large discontinuities in STT_A
% %   %solution: impose a maximum discontinuity of 200km before discard
% %   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   
% %   Flag = [];
% %   for iSeries = 1:1:size(Results.STT_A,1);
% %     dA= diff(find(~isnan(Results.STT_A(iSeries,:))));
% %     if max(dA) > 100; Flag(end+1)= iSeries; end
% %   end
% %   
% %   %identify the date of each of these
% %   Dates = cat(1,Dates,unique(floor(nanmean(Results.Time(Flag,:),2))));

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%
  %problem 1: unphysically large wave amplitudes
  %the number is arbitrsry, but should be sufficient to identify  *class* 
  %of problems rather than the specific cases, so they can be fixed
  %solution: (1) fixed minor bug in cruise-finding routine
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  Flag = [];
  for iSeries = 1:1:size(Results.STT_A,1);
    MaxA = max(Results.STT_A(iSeries,:));
    if MaxA > 30; Flag(end+1)= iSeries; end
  end
  
  %identify the date of each of these
  Dates = cat(1,Dates,unique(floor(nanmean(Results.Time(Flag,:),2))));
   
end

Dates = unique(Dates)';

% % % %%
% % % for iFile = 1:1:14;
% % %   File = [LocalDataDir,'/corwin/IAGOS_st/IAGOS_ST_',num2str(Dates(iFile)),'_v4.mat']
% % %   load(File)
% % %   delete(File)
% % % end