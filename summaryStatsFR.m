function [str1, str2, str3, str4, str5, str6] = summaryStatsFR(area, FR)


%% Fig01
tTest1 = [FR.statsMeanFRPercentile12p5{1}{determineArea(area)}(8).mean1 FR.statsMeanFRPercentile12p5{1}{determineArea(area)}(8).CI1 ...
  FR.statsMeanFRPercentile12p5{1}{determineArea(area)}(8).mean2 FR.statsMeanFRPercentile12p5{1}{determineArea(area)}(8).CI2 ...
  FR.statsMeanFRPercentile12p5{1}{determineArea(area)}(8).fstat FR.statsMeanFRPercentile12p5{1}{determineArea(area)}(8).p];
str1 = sprintf('%s mean(-): %.3g +/- %.3g vs mean(+): %.3g +/- %.3g t = %.3g and p = %.3g \n',...
  area, tTest1(1), tTest1(2), tTest1(3), tTest1(4), tTest1(5), tTest1(6));

wTest1 = [FR.statsVarFRPercentile12p5{1}{determineArea(area)}(8).var1 FR.statsVarFRPercentile12p5{1}{determineArea(area)}(8).var2 ...
  FR.statsVarFRPercentile12p5{1}{determineArea(area)}(8).fstat FR.statsVarFRPercentile12p5{1}{determineArea(area)}(8).p];
str2 = sprintf('%s var(-): %.3g vs var(+): %.3g W = %.3g and p = %.3g \n',...
  area, wTest1(1), wTest1(2), wTest1(3), wTest1(4));


%% Fig02
stats = struct2table(FR.statsVarFRPercentile12p5Split{1}{determineArea(area)}, 'AsArray',true);
iComp1 = stats.iCol1 == 1 & stats.iCol2 == 3;
iComp2 = stats.iCol1 == 2 & stats.iCol2 == 4;
iComp3 = stats.iCol1 == 1 & stats.iCol2 == 2;
iComp4 = stats.iCol1 == 3 & stats.iCol2 == 4;
str3 = ['Constricted(-)(\sigma^2=' num2str(stats.var1(iComp1)) ')vDilated(-)(\sigma^2=' num2str(stats.var2(iComp1)) ') W=' num2str(stats.fstat(iComp1)) ' p=' num2str(stats.p(iComp1)) ', '...
  'Constricted(+)(\sigma^2=' num2str(stats.var1(iComp2)) ')vDilated(+)(\sigma^2=' num2str(stats.var2(iComp2)) ') W=' num2str(stats.fstat(iComp2)) ' p=' num2str(stats.p(iComp2))];
str4 = ['Constricted(-)(\sigma^2=' num2str(stats.var1(iComp3)) ')vConstricted(+)(\sigma^2=' num2str(stats.var2(iComp3)) ') W=' num2str(stats.fstat(iComp3)) ' p=' num2str(stats.p(iComp3)) ', '...
  'Dilated(-)(\sigma^2=' num2str(stats.var1(iComp4)) ')vDilated(+)(\sigma^2=' num2str(stats.var2(iComp4)) ') W=' num2str(stats.fstat(iComp4)) ' p=' num2str(stats.p(iComp4))];

stats = FR.statsMeanFRPercentile12p5Split{1}{determineArea(area)};
str5 = ['Constricted(-)(\mu=' num2str(stats.meansNegative(1)) ' ci=' num2str(stats.ciNegative(1)) ')vDilated(-)(\mu=' num2str(stats.meansNegative(2)) ' ci=' num2str(stats.ciNegative(2)) ') t= ' num2str(stats.statsNegative.tstat) ' p=' num2str(stats.pNegative) ', '...
  'Constricted(+)(\mu=' num2str(stats.meansPositive(1)) ' ci=' num2str(stats.ciPositive(1)) ')vDilated(+)(\mu=' num2str(stats.meansPositive(2)) ' ci=' num2str(stats.ciPositive(2)) ') t= ' num2str(stats.statsPositive.tstat) ' p=' num2str(stats.pPositive)];
str6 = ['Constricted(-)(\mu=' num2str(stats.meansConstricted(1)) ' ci=' num2str(stats.ciConstricted(1)) ')vConstricted(+)(\mu=' num2str(stats.meansConstricted(2)) ' ci=' num2str(stats.ciConstricted(2)) ') t= ' num2str(stats.statsConstricted.tstat) ' p=' num2str(stats.pConstricted) ', '...
  'Dilated(-)(\mu=' num2str(stats.meansDilated(1)) ' ci=' num2str(stats.ciDilated(1)) ')vDilated(+)(\mu=' num2str(stats.meansDilated(2)) ' ci=' num2str(stats.ciDilated(2)) ') t= ' num2str(stats.statsDilated.tstat) ' p=' num2str(stats.pDilated)];