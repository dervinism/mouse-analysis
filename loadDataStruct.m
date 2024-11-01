function dataStruct = loadDataStruct(filePath, type, typeRun)

if nargin < 3
  typeRun = 'run';
end
if nargin < 2
  type = 'all';
end

if strcmpi(typeRun, 'noRun')
  [filePath, fileName, ext] = fileparts(filePath);
  filePath = [filePath '_noRun' filesep fileName ext];
end

dataStruct = load(filePath);
if strcmpi(type,'units')
  if isfield(dataStruct, 'seriesData_ca_muas')
    dataStruct = rmfield(dataStruct, 'seriesData_ca_muas');
  end
  if isfield(dataStruct, 'popData_ca_muas')
    dataStruct = rmfield(dataStruct, 'popData_ca_muas');
  end
  if isfield(dataStruct, 'seriesData_ca_muas_positive')
    dataStruct = rmfield(dataStruct, 'seriesData_ca_muas_positive');
  end
  if isfield(dataStruct, 'popData_ca_muas_positive')
    dataStruct = rmfield(dataStruct, 'popData_ca_muas_positive');
  end
  if isfield(dataStruct, 'seriesData_ca_muas_negative')
    dataStruct = rmfield(dataStruct, 'seriesData_ca_muas_negative');
  end
  if isfield(dataStruct, 'popData_ca_muas_negative')
    dataStruct = rmfield(dataStruct, 'popData_ca_muas_negative');
  end
elseif strcmpi(type,'muas')
  if isfield(dataStruct, 'seriesData_ca')
    dataStruct = rmfield(dataStruct, 'seriesData_ca');
  end
  if isfield(dataStruct, 'popData_ca')
    dataStruct = rmfield(dataStruct, 'popData_ca');
  end
  if isfield(dataStruct, 'seriesData_ca_positive')
    dataStruct = rmfield(dataStruct, 'seriesData_ca_positive');
  end
  if isfield(dataStruct, 'popData_ca_positive')
    dataStruct = rmfield(dataStruct, 'popData_ca_positive');
  end
  if isfield(dataStruct, 'seriesData_ca_negative')
    dataStruct = rmfield(dataStruct, 'seriesData_ca_negative');
  end
  if isfield(dataStruct, 'popData_ca_negative')
    dataStruct = rmfield(dataStruct, 'popData_ca_negative');
  end
end