function fH = phaseGraph(titleStr, FOI, xLims, xTicks, edges, yLims, yTicks, histPhaseOverF, cm, visibility, interpShading, scaleColour)
% Produces a colour-coded phase graph for all FOI. A helper function of
% globalFigs.

if nargin < 12
  scaleColour = false;
end
if nargin < 11
  interpShading = true;
end
if contains(titleStr, 'Vs')
  titleStr = strrep(titleStr, 'Vs', ' wrt ');
end

fH = figure('units', 'normalized', 'position', [0.002, .04, 1/2, .88/1.04], 'Visible', visibility);
histPhaseOverF = [histPhaseOverF zeros(size(histPhaseOverF,1),1)];
histPhaseOverF = [histPhaseOverF; zeros(1,size(histPhaseOverF,2))];
xOriginal = 1:size(histPhaseOverF, 2);
xNew = 1:0.2:size(histPhaseOverF, 2);
FOIInterp = interp1(xOriginal, [FOI FOI(end)+40], xNew);
yOriginal = 1:size(histPhaseOverF, 1);
yNew = 1:0.2:size(histPhaseOverF, 1);
edgesInterp = interp1(yOriginal, edges, yNew);
[X, Y] = meshgrid(xNew, yNew);
if scaleColour
  dd = interp2(histPhaseOverF.^0.5, X, Y);
else
  dd = interp2(histPhaseOverF, X, Y);
end
G = gausswin(5)*gausswin(5)'; G = G/sum(G(:));
h = pcolor(FOIInterp, edgesInterp, conv2(dd, G, 'same'));
h.EdgeColor = 'none';
set(gca,'xscale','log') 
if strcmp('grey', cm)
    colormap(flipud(gray))
elseif strcmp('full', cm)
    colormap(jet)
elseif strcmp('hot', cm)
    colormap(flipud(hot))
elseif strcmp('cold', cm)
    colormap(parula)
end
if interpShading
  shading interp;
end
colorbar;
% ax1 = axesProperties(titleStr, 1, 'normal', 'off', 'w', 'Calibri', 45, 1, 2, [0.01 0.025], 'out', 'off', 'k',...
%   'Frequency (Hz)', xLims, xTicks, 'off', 'k', 'Phase (rad)', yLims, yTicks);
ax1 = axesProperties(titleStr, 1, 'normal', 'off', 'w', 'Calibri', 20, 1, 2, [0.01 0.025], 'out', 'off', 'k',...
  'Frequency (Hz)', xLims, xTicks, 'off', 'k', 'Phase (rad)', yLims, yTicks);
ax1.YTickLabel = {'-2\pi','-3\pi/2','-\pi','-\pi/2','0','\pi/2','\pi','3\pi/2','2\pi'};