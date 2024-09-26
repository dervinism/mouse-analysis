function figH = phaseFigPairMouse(data, eMap, colours, titleStr, figStr, prefix, visibility)
% AnPSD_segs_figs helper function to display phase graphs comparing
% electrode array segment pairs.
% Inputs: data - a structure variable with the following fields:
%           frequencies - a vector of frequencies associated with a pair
%           phase - a vector of segment phase preference relative to
%             another segment in a pair
%           phaseCu - phase upper confidence interval
%           phaseCl - phase lower confidence interval
%           segment1 - electrode number id vector associated with segment 1
%           segment2 - electrode number id vector associated with segment 2
%         eMap - numbered electrode array matrix
%         colours - a matrix of colour code vectors corresponding to segs
%         titleStr - figure title (string)
%         figStr - figure file name (string)
%         prefix - prefix to the figure title and the file name (string)
%         visibility - figure window visble ('on') or not ('off').
%
% Output: figH - figure handle (object).
%

frequencies = data.freq;
phase = bestUnwrap(data.phase);
cu = phase + (data.phase_confU - data.phase);
cl = phase - (data.phase - data.phase_confL);
segments{1} = data.seg1;
segments{2} = data.seg2;

% Draw the phase preference
figH = figure('units', 'normalized', 'position', [0.002, .04, 1, .88], 'Visible', visibility);
semilogx(frequencies, phase, 'r')
hold on
phase(isnan(cu)) = NaN;
semilogx(frequencies, phase, 'k')
semilogx(frequencies, cu, 'k--')
semilogx(frequencies, cl, 'k--')
xlim([frequencies(1) frequencies(end)])
yl = ylim;
ylim([min([-pi yl(1)]) max([pi yl(2)])])

% Estimate dimensions of the electrode array
yl = ylim;
vOffset = 1/36; % proportion of the figure space
height = 34/36;
vPos = [yl(2)-(vOffset)*(yl(2)-yl(1)) yl(2)-(height)*(yl(2)-yl(1))];
xl = log10(xlim);
hOffset = 1/36;
width = 1/27;
hPosLog = [xl(2)-(hOffset)*(xl(2)-xl(1)) xl(2)-(width)*(xl(2)-xl(1))];
hPos = 10.^hPosLog;
arrayXY = [hPos(2) vPos(2) hPos(1)-hPos(2) vPos(1)-vPos(2)];

% Estimate single electrode dimensions
xSide = size(eMap,2);
ySide = size(eMap,1);
exSide = (hPosLog(1)-hPosLog(2))/xSide;
eySide = (vPos(1)-vPos(2))/ySide;

% Draw the electrode array
eGridH = repmat(1:size(eMap,2),size(eMap,1),1);
eGridV = flipud(repmat((1:size(eMap,1))',1,size(eMap,2)));
% for i = 1:size(eGrid,2)
%   addVec = (1:10)' - i;
%   addVec(addVec<0) = 0;
%   eGrid(:,i) = eGrid(:,i) + addVec;
% end
% eGrid = flipud(eGrid);
eMapV = reshape(eMap, 1, size(eMap,1)*size(eMap,2));
eMapV = eMapV(~isnan(eMapV));
for e = 1:numel(eMapV)
  eInd = find(eMap == eMapV(e));
  ePos = [10^(hPosLog(2) + (eGridH(eInd)-1)*exSide) arrayXY(2) + (eGridV(eInd)-1)*eySide...
    10^(hPosLog(2) + (eGridH(eInd)-1)*exSide + exSide) - 10^(hPosLog(2) + (eGridH(eInd)-1)*exSide) eySide];
  rectangle('Position',ePos)
end

% Draw segments
for s = 1:numel(segments)
  segment = segments{s};
  for e = 1:numel(segment)
    eInd = find(eMap == segment(e));
    ePos = [10^(hPosLog(2) + (eGridH(eInd)-1)*exSide) arrayXY(2) + (eGridV(eInd)-1)*eySide...
      10^(hPosLog(2) + (eGridH(eInd)-1)*exSide + exSide) - 10^(hPosLog(2) + (eGridH(eInd)-1)*exSide) eySide];
    rectangle('Position',ePos,'FaceColor',colours(s,:))
  end
end
hold off

title([prefix titleStr], 'Interpreter', 'none')
xlabel('Frequency (Hz)')
ylabel('Phase (Rad)')
prefix = regexprep(prefix,'\','');
saveas(figH, [prefix figStr], 'png')