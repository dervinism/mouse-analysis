function [str1, str2, str3, str4] = summaryStatsFRdots(area, FR)

tickLabels = {'\pi:5\pi/6','5\pi/6:4\pi/6','4\pi/6:3\pi/6','3\pi/6:2\pi/6','2\pi/6:\pi/6','\pi/6:0','0:-\pi/6','-\pi/6:-2\pi/6','-2\pi/6:-3\pi/6','-3\pi/6:-4\pi/6','-4\pi/6:-5\pi/6','-5\pi/6:-\pi'};


%% Fig01
% Var test
data = getLog(FR.areaFRFiltBP0p01to0p05HzIndividual_piOver6{1}{determineArea(area)});
str2 = varsTest(data, tickLabels);

% Means t-test
str1 = meansTest(data, tickLabels);


%% Fig02
% Var test
data = getLog(FR.areaFRFiltBP0p1to0p5HzIndividual_piOver6{1}{determineArea(area)});
str4 = varsTest(data, tickLabels);

% Means t-test
str3 = meansTest(data, tickLabels);



%% Local functions
function str = varsTest(data, tickLabels)

stats = struct2table(varTest(data), 'AsArray',true);
[~, iComp] = max(table2array(stats(:,10)));
dist1 = tickLabels{stats.iCol1(iComp)};
dist2 = tickLabels{stats.iCol2(iComp)};
var1 = round(stats.var1(iComp),2);
var2 = round(stats.var2(iComp),2);
W = table2array(stats(iComp,10));
pval = stats.p(iComp);
str = [dist1 ' (\sigma^2=' num2str(var1) ') vs ' dist2 ' (\sigma^2=' num2str(var2) ') W=' num2str(W) ' p=' num2str(pval)];

function str = meansTest(data, tickLabels)

[dataMean, dataCI95] = datamean(data);
[~, minCol] = min(dataMean);
[~, maxCol] = max(dataMean);
[~, pval, ~, stats] = ttest(data(:,minCol), data(:,maxCol));
dist1 = tickLabels{minCol};
dist2 = tickLabels{maxCol};
mean1 = dataMean(minCol);
mean2 = dataMean(maxCol);
ci1 = dataCI95(2,minCol);
ci2 = dataCI95(2,maxCol);
str = [dist1 ' (\mu=' num2str(mean1) ' ci=' num2str(ci1) ') vs ' dist2 ' (\mu=' num2str(mean2) ' ci=' num2str(ci2) ') t=' num2str(stats.tstat) ' p=' num2str(pval)];