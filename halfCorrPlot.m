function [figPhase, figCoherence, phaseVec, coherenceVec, cohConfVec, phaseCohHalves] = halfCorrPlot(phaseCohHalves, compStr, visibility, varargin)
% Plots correlation graphs for half recording phase and coherence values. A
% helper function for compareHalves_figs, and globalFigs.

if nargin > 3
  uSeries = varargin{1};
else
  uSeries = [];
end
if nargin > 4
  uCentres = varargin{2};
else
  uCentres = [];
end
if nargin > 5
  chColours = varargin{3};
else
  chColours = [];
end
if nargin > 6
  chAreas = varargin{4};
else
  chAreas = [];
end


% INITIALISE FIGURES
if isempty(phaseCohHalves)
  figPhase = [];
  figCoherence = [];
  phaseVec = [];
  coherenceVec = [];
  cohConfVec = [];
  phaseCohHalves = [];
  return
else
  for i = 1:numel(phaseCohHalves)
    if isfield(phaseCohHalves{i}, compStr{1}) && isfield(phaseCohHalves{i}, compStr{2})
      actFOI = phaseCohHalves{i}.(compStr{1}).FOI;
      if ~isempty(actFOI)
        for j = 1:numel(actFOI)
          figPhase(j) = figure('units', 'normalized', 'position', [0.002, .04, 1/2, .88/1.04], 'Visible', visibility); %#ok<*AGROW>
          title(['Phase: ' num2str(actFOI(j)) ' Hz']);
          xlabel([compStr{1} ' phase (rad)'])
          ylabel([compStr{2} ' phase (rad)'])
          figCoherence(j) = figure('units', 'normalized', 'position', [0.002, .04, 1/2, .88/1.04], 'Visible', visibility);
          title(['Coherence: ' num2str(actFOI(j)) ' Hz']);
          xlabel([compStr{1} ' coherence'])
          ylabel([compStr{2} ' coherence'])
        end
        break
      end
    end
  end
end
if ~exist('figPhase','var')
  figPhase = [];
  figCoherence = [];
  phaseVec = [];
  coherenceVec = [];
  cohConfVec = [];
  phaseCohHalves = [];
  return
end


% DRAW FIGURES
xlimMinPhase = zeros(size(figPhase));
xlimMaxPhase = zeros(size(figPhase));
ylimMinPhase = zeros(size(figPhase));
ylimMaxPhase = zeros(size(figPhase));
xlimMinCoherence = zeros(size(figCoherence));
xlimMaxCoherence = zeros(size(figCoherence));
ylimMinCoherence = zeros(size(figCoherence));
ylimMaxCoherence = zeros(size(figCoherence));
phaseVec{j} = [];
coherenceVec{j} = [];
cohConfVec{j} = [];
for i = 1:numel(phaseCohHalves)
  if isfield(phaseCohHalves{i}, compStr{1}) && isfield(phaseCohHalves{i}, compStr{2})
    for j = 1:numel(phaseCohHalves{i}.(compStr{1}).FOI)
      if numel(phaseCohHalves{i}.(compStr{1}).phase) == 1 &&...
          (isnan((phaseCohHalves{i}.(compStr{1}).phase)) || isempty((phaseCohHalves{i}.(compStr{1}).phase)))
        firstValue = nan;
      elseif ~numel(phaseCohHalves{i}.(compStr{1}).phase)
        firstValue = nan;
      else
        firstValue = shiftRadAxis(phaseCohHalves{i}.(compStr{1}).phase(j), [-pi/2 3*pi/2]);
      end
      if numel(phaseCohHalves{i}.(compStr{2}).phase) == 1 &&...
          (isnan((phaseCohHalves{i}.(compStr{2}).phase)) || isempty((phaseCohHalves{i}.(compStr{2}).phase)))
        secondValue = nan;
      elseif ~numel(phaseCohHalves{i}.(compStr{2}).phase)
        secondValue = nan;
      else
        secondValue = shiftRadAxis(phaseCohHalves{i}.(compStr{2}).phase(j), [-pi/2 3*pi/2]);
      end
      figure(figPhase(j));
      hold on
      if ~isnan(firstValue) && ~isnan(secondValue)
        if isempty(uCentres)
          plot(firstValue,secondValue, '.', 'MarkerSize',20)
        else
          for k = 1:numel(chColours{uSeries(i)})
            if sum(uCentres(i) == chColours{uSeries(i)}{k})
              if k == 1
                p1{j} = plot(firstValue,secondValue, 'r.', 'MarkerSize',20);
              elseif k == 2
                p2{j} = plot(firstValue,secondValue, 'b.', 'MarkerSize',20);
              else
                p3{j} = plot(firstValue,secondValue, 'g.', 'MarkerSize',20);
              end
              break
            end
          end
        end
        if xlimMinPhase(j) > min([firstValue secondValue])
          xlimMinPhase(j) = min([firstValue secondValue]);
        end
        if xlimMaxPhase(j) < max([firstValue secondValue])
          xlimMaxPhase(j) = max([firstValue secondValue]);
        end
        xlim([xlimMinPhase(j)-0.05 xlimMaxPhase(j)+0.05]);
        if ylimMinPhase(j) > min([firstValue secondValue])
          ylimMinPhase(j) = min([firstValue secondValue]);
        end
        if ylimMaxPhase(j) < max([firstValue secondValue])
          ylimMaxPhase(j) = max([firstValue secondValue]);
        end
        ylim([ylimMinPhase(j)-0.05 ylimMaxPhase(j)+0.05]);
      end
      phaseVec{j} = [phaseVec{j} [firstValue; secondValue]];
      hold off
      
      figure(figCoherence(j));
      hold on
      if numel(phaseCohHalves{i}.(compStr{1}).coh) == 1 &&...
          (isnan((phaseCohHalves{i}.(compStr{1}).coh)) || isempty((phaseCohHalves{i}.(compStr{1}).coh)))
        firstValue = nan;
        firstValue_conf = nan;
      elseif ~numel(phaseCohHalves{i}.(compStr{1}).coh)
        firstValue = nan;
        firstValue_conf = nan;
      else
        firstValue = phaseCohHalves{i}.(compStr{1}).coh(j);
        firstValue_conf = phaseCohHalves{i}.(compStr{1}).coh_conf(j);
      end
      if numel(phaseCohHalves{i}.(compStr{2}).coh) == 1 &&...
          (isnan((phaseCohHalves{i}.(compStr{2}).coh)) || isempty((phaseCohHalves{i}.(compStr{2}).coh)))
        secondValue = nan;
        secondValue_conf = nan;
      elseif ~numel(phaseCohHalves{i}.(compStr{2}).coh)
        firstValue = nan;
        firstValue_conf = nan;
      else
        secondValue = phaseCohHalves{i}.(compStr{2}).coh(j);
        secondValue_conf = phaseCohHalves{i}.(compStr{2}).coh_conf(j);
      end
      if firstValue > 1
        firstValue = 1;
      end
      if secondValue > 1
        secondValue = 1;
      end
      if ~(firstValue - firstValue_conf <= 0 ||...
          secondValue - secondValue_conf <= 0)
        if ~isnan(firstValue) && ~isnan(secondValue)
          if isempty(uCentres)
            plot(firstValue,secondValue, '.', 'MarkerSize',20)
          else
            for k = 1:numel(chColours{uSeries(i)})
              if sum(uCentres(i) == chColours{uSeries(i)}{k})
                if k == 1
                  p4{j} = plot(firstValue,secondValue, 'r.', 'MarkerSize',20);
                elseif k == 2
                  p5{j} = plot(firstValue,secondValue, 'b.', 'MarkerSize',20);
                else
                  p6{j} = plot(firstValue,secondValue, 'g.', 'MarkerSize',20);
                end
                break
              end
            end
          end
          if xlimMinCoherence(j) > min([firstValue secondValue])
            xlimMinCoherence(j) = min([firstValue secondValue]);
          end
          if xlimMaxCoherence(j) < max([firstValue secondValue])
            xlimMaxCoherence(j) = max([firstValue secondValue]);
          end
          xlim([xlimMinCoherence(j)-0.01 xlimMaxCoherence(j)+0.01]);
          if ylimMinCoherence(j) > min([firstValue secondValue])
            ylimMinCoherence(j) = min([firstValue secondValue]);
          end
          if ylimMaxCoherence(j) < max([firstValue secondValue])
            ylimMaxCoherence(j) = max([firstValue secondValue]);
          end
          ylim([ylimMinCoherence(j)-0.01 ylimMaxCoherence(j)+0.01]);
        end
      end
      coherenceVec{j} = [coherenceVec{j} [firstValue; secondValue]];
      cohConfVec{j} = [cohConfVec{j} [firstValue_conf; secondValue_conf]];
      hold off
    end
  end
end
if exist('p1','var')
  p1{numel(figPhase)+1} = [];
end
if exist('p2','var')
  p2{numel(figPhase)+1} = [];
end
if exist('p3','var')
  p3{numel(figPhase)+1} = [];
end
if exist('p4','var')
  p4{numel(figPhase)+1} = [];
end
if exist('p5','var')
  p5{numel(figPhase)+1} = [];
end
if exist('p6','var')
  p6{numel(figPhase)+1} = [];
end

if exist('p1','var')
  for j = 1:numel(figPhase)
    figure(figPhase(j));
    if (exist('p1','var') && ~isempty(p1{j})) && (exist('p2','var') && ~isempty(p2{j})) && (exist('p3','var') && ~isempty(p3{j}))
      legend([p1{j} p2{j} p3{j}], chAreas{1}, chAreas{2}, chAreas{3}, 'Location','northwest');
    elseif (exist('p1','var') && ~isempty(p1{j})) && (exist('p2','var') && ~isempty(p2{j}))
      legend([p1{j} p2{j}], chAreas{1}, chAreas{2}, 'Location','northwest');
    elseif (exist('p1','var') && ~isempty(p1{j})) && (exist('p3','var') && ~isempty(p3{j}))
      legend([p1{j} p3{j}], chAreas{1}, chAreas{3}, 'Location','northwest');
    elseif (exist('p2','var') && ~isempty(p2{j})) && (exist('p3','var') && ~isempty(p3{j}))
      legend([p2{j} p3{j}], chAreas{2}, chAreas{3}, 'Location','northwest');
    elseif (exist('p1','var') && ~isempty(p1{j}))
      legend(p1{j}, chAreas{1}, 'Location','northwest');
    elseif (exist('p2','var') && ~isempty(p2{j}))
      legend(p2{j}, chAreas{2}, 'Location','northwest');
    elseif (exist('p3','var') && ~isempty(p3{j}))
      legend(p3{j}, chAreas{3}, 'Location','northwest');
    end
    figure(figCoherence(j));
    if (exist('p4','var') && ~isempty(p4{j})) && (exist('p5','var') && ~isempty(p5{j})) && (exist('p6','var') && ~isempty(p6{j}))
      legend([p4{j} p5{j} p6{j}], chAreas{1}, chAreas{2}, chAreas{3}, 'Location','northwest');
    elseif (exist('p4','var') && ~isempty(p4{j})) && (exist('p5','var') && ~isempty(p5{j}))
      legend([p4{j} p5{j}], chAreas{1}, chAreas{2}, 'Location','northwest');
    elseif (exist('p4','var') && ~isempty(p4{j})) && (exist('p6','var') && ~isempty(p6{j}))
      legend([p4{j} p6{j}], chAreas{1}, chAreas{3}, 'Location','northwest');
    elseif (exist('p5','var') && ~isempty(p5{j})) && (exist('p6','var') && ~isempty(p6{j}))
      legend([p5{j} p6{j}], chAreas{2}, chAreas{3}, 'Location','northwest');
    elseif (exist('p4','var') && ~isempty(p4{j}))
      legend(p4{j}, chAreas{1}, 'Location','northwest');
    elseif (exist('p5','var') && ~isempty(p5{j}))
      legend(p5{j}, chAreas{2}, 'Location','northwest');
    elseif (exist('p6','var') && ~isempty(p6{j}))
      legend(p6{j}, chAreas{3}, 'Location','northwest');
    end
  end
end