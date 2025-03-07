% Creates a database for a recording session.

% Set the folder containing all of your animal recording data (top directory)
topDir = uolSourceDirNeuropixels;

% Animal's name or id
animal = 'M191106_MD';

if contains(topDir, 'allensdk')
  repository = 'allensdk';
elseif contains(topDir, 'uol')
  repository = 'uol';
elseif contains(topDir, 'ibl')
  repository = 'ibl';
end

researcher = 'MD';

% Set the configuration
if strcmp(repository, 'allensdk')
  load([topDir filesep animal filesep animal(2:end) '.mat']);
  srRecording = 3e4; % original sampling rate during the recording (approximate only)
  srRecordingLFP = 1250; % LFP sampling rate (approximate only). Applies only to Neuropixels probes.
  [series, ~, nSeries, unitMetadataOI] = determineSeries(headerUnits, unitMetadata); % taylor unit metadata
  shanks = ones(1,nSeries); % number of shanks
  shank = num2cell(ones(nSeries,1)); % shank IDs. For example, if there are 4 shanks, shank = 1:4.
  chN = num2cell(385.*ones(nSeries,1)); % number of recording sites
  LFPch = cellOfMatrices(1:384, nSeries, 1); % LFP channels
  eyeCameraCh = cellOfMatrices([], nSeries, 1); % eye camera channel
  chOI = cellOfMatrices(1:384, nSeries, 1); % channels of interest
  period = cellOfMatrices([0 99999], nSeries, 1); % spiking time period
  periodLFP = cellOfMatrices([0 99999], nSeries, 1); % LFP time period
elseif strcmp(repository, 'uol')
  srRecording = 3e4;
  srRecordingLFP = 2500;
  shanks = ones(1,8);
  shank = {1; 1; 1; 1; 1; 1; 1; 1};
  chN = {385; 385; 385; 385; 385; 385; 385; 385};
  LFPch = {1:384; 1:384; 1:384; 1:384; 1:384; 1:384; 1:384; 1:384};
  eyeCameraCh = {385; 385; 385; 385; 385; 385; 385; 385};
  chOI = {1:30; 31:98; 99:138; 1:138; 139:220; 221:306; 139:306; 307:384};
  period = {[0 99999];[0 99999];[0 99999];[0 99999];[0 99999];[0 99999];[0 99999];[0 99999]};
  periodLFP = period;
elseif strcmp(repository, 'ibl')
  
end

if strcmp(repository, 'allensdk')
  db = allensdkDB(animal(2:end), series);
else
  
  db = [];
  
  i = 1;
  db(i).series = {20191106163201; '201911061632012'};
  db(i).basefilename = 'continuous.imec0.ap_CAR';
  
  i = i+1;
  db(i).series = {20191106163201; '201911061632013'};
  db(i).basefilename = 'continuous.imec0.ap_CAR';
  
  i = i+1;
  db(i).series = {20191106163201; '201911061632014'};
  db(i).basefilename = 'continuous.imec0.ap_CAR';
  
  i = i+1;
  db(i).series = {20191106163201; '2019110616320124'};
  db(i).basefilename = 'continuous.imec0.ap_CAR';
  
  i = i+1;
  db(i).series = {20191106163201; '201911061632015'};
  db(i).basefilename = 'continuous.imec0.ap_CAR';
  
  i = i+1;
  db(i).series = {20191106163201; '201911061632016'};
  db(i).basefilename = 'continuous.imec0.ap_CAR';
  
  i = i+1;
  db(i).series = {20191106163201; '2019110616320156'};
  db(i).basefilename = 'continuous.imec0.ap_CAR';
  
  i = i+1;
  db(i).series = {20191106163201; '201911061632017'};
  db(i).basefilename = 'continuous.imec0.ap_CAR';
end


db = commonDB(db, animal, repository, researcher, shanks, shank,...
  chN, LFPch, eyeCameraCh, chOI, period, periodLFP);


for i = 1:length(db)
  if ~isfield(db(i), 'series') || isempty(db(i).series)
    db(i).series = 1;
  end
  if ~isfield(db(i), 'entryName') || isempty(db(i).entryName)
    if iscell(db(i).series) && ~ischar(db(i).series{2})
      seriesFirst = num2str(db(i).series{1});
      seriesLast = num2str(db(i).series{end});
      db(i).entryName = [db(i).animal '_s' seriesFirst seriesLast(end)];
    elseif iscell(db(i).series) && ischar(db(i).series{2})
      db(i).entryName = [db(i).animal '_s' db(i).series{2}];
    else
      db(i).entryName = [db(i).animal '_s' num2str(db(i).series)];
    end
  end
end


clear animal chN eyeCameraCh i LFPch chOI shank shanks nSeries period periodLFP researcher
fprintf('Created db\n')
