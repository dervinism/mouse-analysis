function [s, sa] = combineMotionData(ifilename, ofilename)
% Combines two or more total movement data analysis files.

for iFile = 1:numel(ifilename)
  data{iFile} = load(ifilename{iFile});
end

s = data{1}.s;
sa = data{1}.sa;
for iFile = 2:numel(ifilename)
  s = [s data{iFile}.s]; %#ok<*AGROW>
  sa = [sa data{iFile}.sa];
end

save(ofilename, 's', 'sa');