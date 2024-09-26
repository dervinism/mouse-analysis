function figH = cohMatrix(matFile, iF, condition12, area1, area2, condition34, area3, area4) %#ok<*STOUT>
% figH = cohMatrix(matFile, iF, condition, area1, area2, area3, area4)
% Function for displaying coherence matrix, other related matrices, and
% overlaps. Function performs overlap calculation and paired-sample t-test
% only if the the two comparisons are coming from the same condition.
% Inputs: matFile - mat file name with pre-processed data (area2areaCohMats.mat);
%         iF - frequency of interest;
%         condition12 - 'awake', 'anaesthesia', or 'all' fo the first comparison corresponding to areas 1 and 2;
%         area1, area2 - the first comparison: area1 vs area2 (for overlap calculation);
%         condition34 - 'awake', 'anaesthesia', or 'all' fo the second comparison corresponding to areas 3 and 4;
%         area3, area4 - the second comparison: area3 vs area4 (for overlap calculation).
% Output: figH - figure handle to coherence colour map (not set up yet).

load(matFile); %#ok<*LOAD>

% PERFORM INPUT CHECKS
indFOI = find(FOI == iF);
if isempty(indFOI)
  error('No values exist for your supplied frequency of interest');
end

for iCond = 1:numel(conditions) %#ok<*USENS>
  if strcmpi(conditions{iCond}, condition12) %#ok<*IDISVAR>
    conditionInd12 = iCond;
  end
end
if ~exist('conditionInd12','var')
  error('No values exist for your supplied condition');
end

for iCond = 1:numel(conditions) %#ok<*USENS>
  if strcmpi(conditions{iCond}, condition34) %#ok<*IDISVAR>
    conditionInd34 = iCond;
  end
end
if ~exist('conditionInd34','var')
  error('No values exist for your supplied condition');
end

% CONCATENATE COHERENCE AND OTHER MATRICES FOR YOUR FOI
recMat_cond_iF12 = NaN(numel(comparisonAreas),numel(comparisonAreas),numel(cohMat{conditionInd12}));
cohMat_cond_iF12 = NaN(numel(comparisonAreas),numel(comparisonAreas),numel(cohMat{conditionInd12}));
mfrPR1Mat_cond_iF12 = NaN(numel(comparisonAreas),numel(comparisonAreas),numel(cohMat{conditionInd12}));
mfrPR2Mat_cond_iF12 = NaN(numel(comparisonAreas),numel(comparisonAreas),numel(cohMat{conditionInd12}));
signMat_cond_iF12 = NaN(numel(comparisonAreas),numel(comparisonAreas),numel(cohMat{conditionInd12}));
numberOfRecordings12 = numel(cohMat{conditionInd12});
for iRec = 1:numberOfRecordings12
  recMat_cond_iF12(:,:,iRec) = recMat{conditionInd12}{iRec}(:,:,indFOI);
  cohMat_cond_iF12(:,:,iRec) = cohMat{conditionInd12}{iRec}(:,:,indFOI);
  mfrPR1Mat_cond_iF12(:,:,iRec) = mfrPR1Mat{conditionInd12}{iRec}(:,:,indFOI);
  mfrPR2Mat_cond_iF12(:,:,iRec) = mfrPR2Mat{conditionInd12}{iRec}(:,:,indFOI);
  signMat_cond_iF12(:,:,iRec) = ~isnan(cohMat{conditionInd12}{iRec}(:,:,indFOI));
end

recMat_cond_iF34 = NaN(numel(comparisonAreas),numel(comparisonAreas),numel(recMat{conditionInd34}));
cohMat_cond_iF34 = NaN(numel(comparisonAreas),numel(comparisonAreas),numel(cohMat{conditionInd34}));
mfrPR1Mat_cond_iF34 = NaN(numel(comparisonAreas),numel(comparisonAreas),numel(cohMat{conditionInd34}));
mfrPR2Mat_cond_iF34 = NaN(numel(comparisonAreas),numel(comparisonAreas),numel(cohMat{conditionInd34}));
signMat_cond_iF34 = NaN(numel(comparisonAreas),numel(comparisonAreas),numel(cohMat{conditionInd34}));
numberOfRecordings34 = numel(cohMat{conditionInd34});
for iRec = 1:numberOfRecordings34
  recMat_cond_iF34(:,:,iRec) = recMat{conditionInd34}{iRec}(:,:,indFOI);
  cohMat_cond_iF34(:,:,iRec) = cohMat{conditionInd34}{iRec}(:,:,indFOI);
  mfrPR1Mat_cond_iF34(:,:,iRec) = mfrPR1Mat{conditionInd34}{iRec}(:,:,indFOI);
  mfrPR2Mat_cond_iF34(:,:,iRec) = mfrPR2Mat{conditionInd34}{iRec}(:,:,indFOI);
  signMat_cond_iF34(:,:,iRec) = ~isnan(cohMat{conditionInd34}{iRec}(:,:,indFOI));
end

% CALCULATE OVERLAP BETWEEN ALL COMPARISON CONDITIONS
combos = combntns(1:numel(comparisonAreas),2);
combosOfCombos = combntns(1:size(combos,1),2);

overlapVec12 = zeros(size(combosOfCombos,1),numberOfRecordings12);
overlapSum12 = zeros(size(combosOfCombos,1),1);
cohMat_cond_iF12_logical = ~isnan(cohMat_cond_iF12);
for iComb = 1:size(combosOfCombos,1)
  vec1 = cohMat_cond_iF12_logical(combos(combosOfCombos(iComb,1),1), combos(combosOfCombos(iComb,1),2), :);
  vec1 = reshape(vec1,1,numel(vec1));
  vec2 = cohMat_cond_iF12_logical(combos(combosOfCombos(iComb,2),1), combos(combosOfCombos(iComb,2),2), :);
  vec2 = reshape(vec2,1,numel(vec2));
  overlapVec12(iComb,:) = vec1 & vec2;
  overlapSum12(iComb) = sum(vec1 & vec2);
end

overlapVec34 = zeros(size(combosOfCombos,1),numberOfRecordings34);
overlapSum34 = zeros(size(combosOfCombos,1),1);
cohMat_cond_iF34_logical = ~isnan(cohMat_cond_iF34);
for iComb = 1:size(combosOfCombos,1)
  vec1 = cohMat_cond_iF34_logical(combos(combosOfCombos(iComb,1),1), combos(combosOfCombos(iComb,1),2), :);
  vec1 = reshape(vec1,1,numel(vec1));
  vec2 = cohMat_cond_iF34_logical(combos(combosOfCombos(iComb,2),1), combos(combosOfCombos(iComb,2),2), :);
  vec2 = reshape(vec2,1,numel(vec2));
  overlapVec34(iComb,:) = vec1 & vec2;
  overlapSum34(iComb) = sum(vec1 & vec2);
end

% DRAW A MATRIX FOR YOUR FOI
% figH = figure;
% h = pcolor(1:numel(comparisonAreas), 1:numel(comparisonAreas), cohMatMean{conditionInd}(:,:,indFOI));
% h.EdgeColor = 'none';
% colormap(parula)
% colorbar;
% ax1 = axesProperties(titleStr, 1, 'normal', 'off', 'w', 'Calibri', 20, 1, 2, [0.01 0.025], 'out', 'off', 'k',...
%   'Frequency (Hz)', xLims, xTicks, 'off', 'k', 'Phase (rad)', yLims, yTicks);
% ax1.YTickLabel = {'-\pi/2','0','\pi/2','\pi','-\pi/2'};

% PRINT OUT THE MATRICES
% Matrices for the first comparison
disp(['Coherence matrix for the first comparison for frequency: ' num2str(iF)]);
disp(cohMatMean{conditionInd12}(:,:,indFOI));
disp('Significant recordings in each matrix entry');
signMat12 = tallyMat{conditionInd12}(:,:,indFOI); %#ok<*NODEF>
% for iDiag = 1:numel(comparisonAreas)
%   signMat(iDiag,iDiag) = 0;
% end
disp(signMat12);

numberOfRecordingsIndividually12 = zeros(numel(comparisonAreas),numel(comparisonAreas));
for iArea = 1:numel(areaRecCount{conditionInd12})
  [iRow, iCol] = areaEntryFromArea(iArea);
  numberOfRecordingsIndividually12(iRow, iCol) = areaRecCount{conditionInd12}{iArea};
  numberOfRecordingsIndividually12(iCol, iRow) = areaRecCount{conditionInd12}{iArea};
end
disp('Total number of recordings in each matrix entry');
disp(numberOfRecordingsIndividually12);

disp(['Total number of recordings: ' num2str(numberOfRecordings12)]);

% Matrices for the second comparison
disp(['Coherence matrix for the second comparison for frequency: ' num2str(iF)]);
disp(cohMatMean{conditionInd34}(:,:,indFOI));
disp('Significant recordings in each matrix entry');
signMat34 = tallyMat{conditionInd34}(:,:,indFOI); %#ok<*NODEF>
% for iDiag = 1:numel(comparisonAreas)
%   signMat(iDiag,iDiag) = 0;
% end
disp(signMat34);

numberOfRecordingsIndividually34 = zeros(numel(comparisonAreas),numel(comparisonAreas));
for iArea = 1:numel(areaRecCount{conditionInd34})
  [iRow, iCol] = areaEntryFromArea(iArea);
  numberOfRecordingsIndividually34(iRow, iCol) = areaRecCount{conditionInd34}{iArea};
  numberOfRecordingsIndividually34(iCol, iRow) = areaRecCount{conditionInd34}{iArea};
end
disp('Total number of recordings in each matrix entry');
disp(numberOfRecordingsIndividually34);

disp(['Total number of recordings: ' num2str(numberOfRecordings34)]);

% Overlaps
if strcmpi(condition12, condition34)
  [overlapCOI, overlapCOIvec] = findOverlap(overlapSum12, overlapVec12, combos, combosOfCombos, area1, area2, area3, area4);
  fprintf('Overlaping significant recordings between %svs%s and %svs%s: %g\n', area1, area2, area3, area4, overlapCOI);
  disp('Recordings that overlap: ');
  disp(overlapCOIvec);
end

% Firing rates
entry1 = areaEntry(area1);
entry2 = areaEntry(area2);
disp('Significant recordings in the first comparison:');
disp(reshape(signMat_cond_iF12(entry2,entry1,:),1,numberOfRecordings12));
disp('Firing rates for area1 in the first comparison:');
disp(reshape(mfrPR1Mat_cond_iF12(entry2,entry1,:),1,numberOfRecordings12));
disp('Firing rates for area2 in the first comparison:');
disp(reshape(mfrPR2Mat_cond_iF12(entry2,entry1,:),1,numberOfRecordings12));

entry3 = areaEntry(area3);
entry4 = areaEntry(area4);
disp('Significant recordings in the second comparison:');
disp(reshape(signMat_cond_iF34(entry4,entry3,:),1,numberOfRecordings34));
disp('Firing rates for area1 in the second comparison:');
disp(reshape(mfrPR1Mat_cond_iF34(entry4,entry3,:),1,numberOfRecordings34));
disp('Firing rates for area2 in the second comparison:');
disp(reshape(mfrPR2Mat_cond_iF34(entry4,entry3,:),1,numberOfRecordings34));

% statistical tests
% Independent-samples t-test for coherence
if strcmpi(condition12, condition34)
  [~, pval] = ttest2(reshape(cohMat_cond_iF12(entry2,entry1,:),1,numberOfRecordings12), reshape(cohMat_cond_iF34(entry4,entry3,:),1,numberOfRecordings34));
  fprintf('p-value of the t-test for coherence comparisons %svs%s and %svs%s: %g\n', area1, area2, area3, area4, pval);
end

% Chi-square test for significant recordings
[~, pval] = chi2generate(signMat12(entry2,entry1), numberOfRecordingsIndividually12(entry2,entry1),...
  signMat34(entry4,entry3), numberOfRecordingsIndividually34(entry4,entry3));
fprintf('p-value of chi2 test for significant recordings comparisons %svs%s and %svs%s: %g\n', area1, area2, area3, area4, pval);

% Independent-samples t-test for firing rates comparing significant and non-significant recordings
recMat_cond_iF12(isnan(recMat_cond_iF12)) = 0;
mfrPRMat_cond_iF12_min = NaN(size(mfrPR1Mat_cond_iF12));
for dim1 = 1:size(mfrPRMat_cond_iF12_min,1)
  for dim2 = 1:size(mfrPRMat_cond_iF12_min,2)
    for dim3 = 1:size(mfrPRMat_cond_iF12_min,3)
      if dim1 > dim2
        mfrPRMat_cond_iF12_min(dim1,dim2,dim3) = min(mfrPR1Mat_cond_iF12(dim1,dim2,dim3), mfrPR2Mat_cond_iF12(dim1,dim2,dim3));
      else
        signMat_cond_iF12(dim1,dim2,dim3) = 0;
        recMat_cond_iF12(dim1,dim2,dim3) = 0;
      end
    end
  end
end
mfrSignRec_min = mfrPRMat_cond_iF12_min(logical(signMat_cond_iF12));
nonsignMat_cond_iF12 = recMat_cond_iF12;
nonsignMat_cond_iF12(logical(signMat_cond_iF12)) = 0;
mfrNonsignRec_min = mfrPRMat_cond_iF12_min(logical(nonsignMat_cond_iF12));
mfrRec_min12 = [mfrSignRec_min; mfrNonsignRec_min];
[~, pval] = ttest2(mfrSignRec_min, mfrNonsignRec_min);
fprintf('p-value of the t-test for firing rate comparisons between significant (MFR %g) and non-significant (MFR %g) recordings for the 1st condition: %g\n',...
  mean(mfrSignRec_min), mean(mfrNonsignRec_min), pval);

recMat_cond_iF34(isnan(recMat_cond_iF34)) = 0;
mfrPRMat_cond_iF34_min = NaN(size(mfrPR1Mat_cond_iF34));
for dim1 = 1:size(mfrPRMat_cond_iF34_min,1)
  for dim2 = 1:size(mfrPRMat_cond_iF34_min,2)
    for dim3 = 1:size(mfrPRMat_cond_iF34_min,3)
      if dim1 > dim2
        mfrPRMat_cond_iF34_min(dim1,dim2,dim3) = min(mfrPR1Mat_cond_iF34(dim1,dim2,dim3), mfrPR2Mat_cond_iF34(dim1,dim2,dim3));
      else
        signMat_cond_iF34(dim1,dim2,dim3) = 0;
        recMat_cond_iF34(dim1,dim2,dim3) = 0;
      end
    end
  end
end
mfrSignRec_min = mfrPRMat_cond_iF34_min(logical(signMat_cond_iF34));
nonsignMat_cond_iF34 = recMat_cond_iF34;
nonsignMat_cond_iF34(logical(signMat_cond_iF34)) = 0;
mfrNonsignRec_min = mfrPRMat_cond_iF34_min(logical(nonsignMat_cond_iF34));
mfrRec_min34 = [mfrSignRec_min; mfrNonsignRec_min];
[~, pval] = ttest2(mfrSignRec_min, mfrNonsignRec_min);
fprintf('p-value of the t-test for firing rate comparisons between significant (MFR %g) and non-significant (MFR %g) recordings for the 2nd condition: %g\n',...
  mean(mfrSignRec_min), mean(mfrNonsignRec_min), pval);

% Independent-samples t-test for firing rates comparing the two conditions
[~, pval] = ttest2(mfrRec_min12, mfrRec_min34);
fprintf('p-value of the t-test for firing rate comparisons between condition 1 (MFR %g) and condition 2 (MFR %g): %g\n',...
  mean(mfrRec_min12), mean(mfrRec_min34), pval);
end


function [overlapCOI, overlapCOIvec] = findOverlap(overlapSum, overlapVec, combos, combosOfCombos, area1, area2, area3, area4)

entry1 = areaEntry(area1);
entry2 = areaEntry(area2);
iCombos1 = findCombo(entry1, entry2, combos);
entry3 = areaEntry(area3);
entry4 = areaEntry(area4);
iCombos2 = findCombo(entry3, entry4, combos);
iCombo = findCombo(iCombos1, iCombos2, combosOfCombos);
overlapCOI = overlapSum(iCombo);
overlapCOIvec = overlapVec(iCombo,:);
end

function iCombos1 = findCombo(entry1, entry2, combos)

if entry1 < entry2
  combos1 = combos(:,1) == entry1;
  combos2 = combos(:,2) == entry2;
elseif entry1 > entry2
  combos1 = combos(:,1) == entry2;
  combos2 = combos(:,2) == entry1;
end
iCombos1 = find(combos1 & combos2);
end

function [iRow, iCol] = areaEntryFromArea(area)

if area == 1
  iRow = 1; iCol = 2; 
elseif area == 2
  iRow = 1; iCol = 5; 
elseif area == 3
  iRow = 1; iCol = 6; 
elseif area == 4
  iRow = 1; iCol = 10; 
elseif area == 5
  iRow = 1; iCol = 7; 
elseif area == 6
  iRow = 1; iCol = 9; 
elseif area == 7
  iRow = 1; iCol = 8;
elseif area == 8
  iRow = 2; iCol = 5; 
elseif area == 9
  iRow = 2; iCol = 6; 
elseif area == 10
  iRow = 2; iCol = 10; 
elseif area == 11
  iRow = 2; iCol = 7; 
elseif area == 12
  iRow = 2; iCol = 9; 
elseif area == 13
  iRow = 2; iCol = 8; 
elseif area == 14
  iRow = 3; iCol = 4; 
elseif area == 15
  iRow = 3; iCol = 5; 
elseif area == 16
  iRow = 3; iCol = 6; 
elseif area == 17
  iRow = 3; iCol = 10; 
elseif area == 18
  iRow = 3; iCol = 7; 
elseif area == 19
  iRow = 3; iCol = 9; 
elseif area == 20
  iRow = 3; iCol = 8;
elseif area == 21
  iRow = 4; iCol = 5; 
elseif area == 22
  iRow = 4; iCol = 6; 
elseif area == 23
  iRow = 4; iCol = 10; 
elseif area == 24
  iRow = 4; iCol = 7; 
elseif area == 25
  iRow = 4; iCol = 9; 
elseif area == 26
  iRow = 4; iCol = 8;
elseif area == 27
  iRow = 5; iCol = 6; 
elseif area == 28
  iRow = 5; iCol = 10; 
elseif area == 29
  iRow = 5; iCol = 7; 
elseif area == 30
  iRow = 5; iCol = 9; 
elseif area == 31
  iRow = 5; iCol = 8;
elseif area == 32
  iRow = 6; iCol = 10; 
elseif area == 33
  iRow = 6; iCol = 7; 
elseif area == 34
  iRow = 6; iCol = 9; 
elseif area == 35
  iRow = 6; iCol = 8;
elseif area == 36
  iRow = 9; iCol = 8;
elseif area == 37
  iRow = 9; iCol = 7;
elseif area == 38
  iRow = 8; iCol = 7;
elseif area == 39
  iRow = 11; iCol = 5; 
elseif area == 40
  iRow = 11; iCol = 6; 
elseif area == 41
  iRow = 11; iCol = 10; 
elseif area == 42
  iRow = 11; iCol = 7; 
elseif area == 43
  iRow = 11; iCol = 9; 
elseif area == 44
  iRow = 11; iCol = 8;
elseif area == 45
  iRow = 12; iCol = 5; 
elseif area == 46
  iRow = 12; iCol = 6; 
elseif area == 47
  iRow = 12; iCol = 10; 
elseif area == 48
  iRow = 12; iCol = 7; 
elseif area == 49
  iRow = 12; iCol = 9; 
elseif area == 50
  iRow = 12; iCol = 8;
elseif area == 51
  iRow = 11; iCol = 13;
elseif area == 52
  iRow = 12; iCol = 13;
elseif area == 53
  iRow = 5; iCol = 13;
end
end

function [chi2stat, pval] = chi2generate(n1, N1, n2, N2)

x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)]; %#ok<*REPMAT>
[~,chi2stat,pval] = crosstab(x1,x2);
end