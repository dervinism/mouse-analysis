options.saturationMethod = 'combined';
options.saturationPlot = true;
options.powerCalcMethod = 'filter';

file = 'R:\Neuropix\Shared\Data\M191106_MD\20191106163201\continuous.imec0.lf.bin';
nCh = 385;

options.chOI = round(nCh/33):round(nCh/33):nCh;
power = lfpPowersPlot(file, nCh, 2500, options);



filePR = 'R:\CSN\md406\runNSG\runNSG_M191106_MD\M191106_MD.mat';
load(filePR);
[~, series] = fileparts(fileparts(file));
fnsData = fieldnames(dataStruct.seriesData);
seriesOI = find(contains(fnsData, series));
for iSeries = 1:numel(seriesOI)
  PR = sum(dataStruct.seriesData.(fnsData{seriesOI(iSeries)}).popData.MUAsAll,1);
  ch = dataStruct.seriesData.(fnsData{seriesOI(iSeries)}).db(seriesOI(iSeries)).chOI;
  time = (1:numel(PR)).*(1/400);
  figure; plot(time,PR);
  title([fnsData{seriesOI(iSeries)} ' ch' num2str(ch(1)) '-' num2str(ch(end))], 'Interpreter','None');
end