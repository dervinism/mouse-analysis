options.saturationMethod = 'combined';
options.saturationPlot = true;


file1 = 'R:\Neuropix\Shared\Data\allensdk\M766640955\probe_773592315_lfp.csv';
nCh1 = 88;

file2 = 'R:\Neuropix\Shared\Data\allensdk\M766640955\probe_773592318_lfp.csv';
nCh2 = 89;

file3 = 'R:\Neuropix\Shared\Data\allensdk\M766640955\probe_773592320_lfp.csv';
nCh3 = 82;

file4 = 'R:\Neuropix\Shared\Data\allensdk\M766640955\probe_773592324_lfp.csv';
nCh4 = 76;

file5 = 'R:\Neuropix\Shared\Data\allensdk\M766640955\probe_773592328_lfp.csv';
nCh5 = 88;

file6 = 'R:\Neuropix\Shared\Data\allensdk\M766640955\probe_773592330_lfp.csv';
nCh6 = 94;

file = file6;
nCh = nCh6;
options.chOI = round(nCh/10):round(nCh/10):nCh;
lfpPowersPlot(file, nCh, 1250, options);



% [~, ~, chIDs1] = loadcsv(file1);
% [~, ~, chIDs2] = loadcsv(file2);
% [~, ~, chIDs3] = loadcsv(file3);
% [~, ~, chIDs4] = loadcsv(file4);
% [~, ~, chIDs5] = loadcsv(file5);
% [~, ~, chIDs6] = loadcsv(file6);
% chIDs = [chIDs1; chIDs2; chIDs3; chIDs4; chIDs5; chIDs6];
% 
% load('R:\Neuropix\Shared\Data\allensdk\M766640955\766640955.mat');
% chIDs_spikes = unique(cell2mat(unitMetadata(:,14)));
% sum(ismember(chIDs_spikes, chIDs))
% numel(chIDs_spikes)