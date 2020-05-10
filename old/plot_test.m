subplot(2,1,1)


load out; a = nansum(Results.Sum,3); b = nansum(Results.Count,3); c = a./b; 

% load out; c = log10(nansum(Results.Count,3));

% load out; i = 9;a = Results.Sum(:,:,i); b = Results.Count(:,:,i); c = a./b; 

c = log10(c);

c = smoothn(c,[3,3]);

m_proj('lambert','lon',[-125,2],'lat',[20,80]);
m_pcolor(xi,yi,c);

m_grid;
m_coast('color','k');
colorbar
caxis([-2.5 -1.7])
% caxis([0 0.02])
redyellowblue16


subplot(2,1,2)


% load out; a = nansum(Results.Sum,3); b = nansum(Results.Count,3); c = a./b; 

load out; c = log10(nansum(Results.Count,3));

% load out; i = 9;a = Results.Sum(:,:,i); b = Results.Count(:,:,i); c = a./b; 

% c = log10(c);

% c = smoothn(c,[3,3]);

m_proj('Lambert','lon',[-125,2],'lat',[20,80]);
m_pcolor(xi,yi,c);

m_grid;
m_coast('color','k');
colorbar
% caxis([-2.5 -1.7])
% caxis([0 0.02])
redyellowblue16