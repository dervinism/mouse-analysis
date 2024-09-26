function entry = areaEntry(areaStr)
% The function provides the area entry number given the area ID. It is a
% helper function of globalPhaseAnalysis and cohMatrix.

if strcmpi(areaStr, '2')       % VB1
  entry = 1;
elseif strcmpi(areaStr, '8')   % VB2
  entry = 2;
elseif strcmpi(areaStr, '24')  % Th1
  entry = 3;
elseif strcmpi(areaStr, '810') % Th2
  entry = 4;
elseif strcmpi(areaStr, '1')   % S1
  entry = 5;
elseif strcmpi(areaStr, '7')   % RSC
  entry = 6;
elseif strcmpi(areaStr, '5')   % DG
  entry = 7;
elseif strcmpi(areaStr, '11')  % CA3
  entry = 8;
elseif strcmpi(areaStr, '6')   % CA1
  entry = 9;
elseif strcmpi(areaStr, '56')  % Hp
  entry = 10;
end
end