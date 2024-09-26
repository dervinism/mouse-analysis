function [figMFR, mfrVec] = halfCorrPlotMFR(mfrHalves, compStr, visibility, uSeries, uCentres, chColours, chAreas)
% Plots correlation graphs for half recording mean firing rate values. A
% helper function for globalFigs.


% INITIALISE FIGURES
figMFR = figure('units', 'normalized', 'position', [0.002, .04, 1/2, .88/1.04], 'Visible', visibility); %#ok<*AGROW>
title('MFR');
xlabel([compStr{1} ' mean firing rate (APs/s)'])
ylabel([compStr{2} ' mean firing rate (APs/s)'])


% DRAW FIGURES
xlimMinMFR = zeros(size(figMFR));
xlimMaxMFR = zeros(size(figMFR));
ylimMinMFR = zeros(size(figMFR));
ylimMaxMFR = zeros(size(figMFR));
mfrVec = [];
for i = 1:numel(mfrHalves)
  if isfield(mfrHalves{i}, compStr{1}) && isfield(mfrHalves{i}, compStr{2})
    fig = figure(figMFR);
    hold on
    if ~isnan(mfrHalves{i}.(compStr{1})) && ~isnan(mfrHalves{i}.(compStr{2}))
      if isempty(uCentres)
        plot(mfrHalves{i}.(compStr{1}),mfrHalves{i}.(compStr{2}), '.', 'MarkerSize',20)
      else
        for k = 1:numel(chColours{uSeries(i)})
          if sum(uCentres(i) == chColours{uSeries(i)}{k})
            if k == 1
              p1 = plot(mfrHalves{i}.(compStr{1}),mfrHalves{i}.(compStr{2}), 'r.', 'MarkerSize',20);
            elseif k == 2
              p2 = plot(mfrHalves{i}.(compStr{1}),mfrHalves{i}.(compStr{2}), 'b.', 'MarkerSize',20);
            else
              p3 = plot(mfrHalves{i}.(compStr{1}),mfrHalves{i}.(compStr{2}), 'g.', 'MarkerSize',20);
            end
            break
          end
        end
      end
      if xlimMinMFR > min([mfrHalves{i}.(compStr{1}) mfrHalves{i}.(compStr{2})])
        xlimMinMFR = min([mfrHalves{i}.(compStr{1}) mfrHalves{i}.(compStr{2})]);
      end
      if xlimMaxMFR < max([mfrHalves{i}.(compStr{1}) mfrHalves{i}.(compStr{2})])
        xlimMaxMFR = max([mfrHalves{i}.(compStr{1}) mfrHalves{i}.(compStr{2})]);
      end
      xlim([xlimMinMFR-0.05 xlimMaxMFR+0.05]);
      if ylimMinMFR > min([mfrHalves{i}.(compStr{1}) mfrHalves{i}.(compStr{2})])
        ylimMinMFR = min([mfrHalves{i}.(compStr{1}) mfrHalves{i}.(compStr{2})]);
      end
      if ylimMaxMFR < max([mfrHalves{i}.(compStr{1}) mfrHalves{i}.(compStr{2})])
        ylimMaxMFR = max([mfrHalves{i}.(compStr{1}) mfrHalves{i}.(compStr{2})]);
      end
      ylim([ylimMinMFR-0.05 ylimMaxMFR+0.05]);
    end
    mfrVec = [mfrVec [mfrHalves{i}.(compStr{1}); mfrHalves{i}.(compStr{2})]];
    hold off
  end
end

if ~isempty(chAreas)
  figure(fig);
  if (exist('p1','var') && ~isempty(p1)) && (exist('p2','var') && ~isempty(p2)) && (exist('p3','var') && ~isempty(p3))
    legend([p1 p2 p3], chAreas{1}, chAreas{2}, chAreas{3}, 'Location','northwest');
  elseif (exist('p1','var') && ~isempty(p1)) && (exist('p2','var') && ~isempty(p2))
    legend([p1 p2], chAreas{1}, chAreas{2}, 'Location','northwest');
  elseif (exist('p1','var') && ~isempty(p1)) && (exist('p3','var') && ~isempty(p3))
    legend([p1 p3], chAreas{1}, chAreas{3}, 'Location','northwest');
  elseif (exist('p2','var') && ~isempty(p2)) && (exist('p3','var') && ~isempty(p3))
    legend([p2 p3], chAreas{2}, chAreas{3}, 'Location','northwest');
  elseif (exist('p1','var') && ~isempty(p1))
    legend(p1, chAreas{1}, 'Location','northwest');
  elseif (exist('p2','var') && ~isempty(p2))
    legend(p2, chAreas{2}, 'Location','northwest');
  elseif (exist('p3','var') && ~isempty(p3))
    legend(p3, chAreas{3}, 'Location','northwest');
  end
end