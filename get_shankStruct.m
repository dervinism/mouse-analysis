function [shankStruct, shank, units, unitMetadata, xcoords, ycoords, spk,...
  MUAs, phaseCoh] = get_shankStruct(dbStruct, sh)

shankIDs = fieldnames(dbStruct.shankData);
shankStruct = eval(['dbStruct.shankData.' shankIDs{sh}]);
if isempty(shankStruct)
  shankStruct = [];
  shank = [];
  units = [];
  unitMetadata = [];
  xcoords = [];
  ycoords = [];
  spk = [];
  MUAs = [];
  phaseCoh = [];
  return
end

shank = shankStruct.shankID;
units = shankStruct.units;
if isempty(units)
  unitMetadata = []; xcoords = []; ycoords = [];
  spk = []; MUAs = [];
  phaseCoh = [];
  return
end
unitMetadata = shankStruct.unitMetadata;
xcoords = shankStruct.unitMetadata(:,4);
ycoords = shankStruct.unitMetadata(:,5);

if isfield(shankStruct, 'spk')
  spk = shankStruct.spk;
else
  spk = [];
end
if isfield(shankStruct, 'MUAs')
  MUAs = shankStruct.MUAs;
else
  MUAs = [];
end

if isfield(shankStruct, 'phaseCoh')
  phaseCoh = shankStruct.phaseCoh;
else
  phaseCoh = [];
end