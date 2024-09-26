% Creates a database for a recording session.
% No inputs. Outputs (self-explanatory): topDir, db, srRecording.

% Set the folder containing all of your animal recording data (top directory)
if ispc && ~isempty(dir('R:\'))
  topDir = 'R:\Neuropix\Shared\Data';
elseif ispc && isempty(dir('R:\'))
  error('Could not find R: drive on your system')
else 
  error('Make sure to set up topDir containing recording data')
end

% Animal's name or id
animal = 'M191106_MD';

% Set the configuration
srRecording = 3e4;                                      % original sampling rate during the recording
srRecordingLFP = 2500;                                  % LFP sampling rate. Applies only to Neuropixels probes.
shanks = ones(1,8);                                     % number of shanks
shank = {1; 1; 1; 1; 1; 1; 1; 1};                       % shank IDs. For example, if there are 4 shanks, shank = 1:4.
V1shank = {0; 0; 0; 0; 0; 0; 0; 0};                     % Is any of the shanks in V1?
shankGroups = {{1}; {1}; {1}; {1}; {1}; {1}; {1}; {1}}; % shank ID groups. For example, if there are 2 shanks in each of two groups, shankGroups = {1:2, 3:4}.
chN = {385; 385; 385; 385; 385; 385; 385; 385};         % number of recording sites
LFPch = {1:384; 1:384; 1:384; 1:384; 1:384; 1:384; 1:384; 1:384};
eyeCameraCh = {385; 385; 385; 385; 385; 385; 385; 385};
chOI = {1:30; 31:98; 99:138; 1:138; 139:220; 221:306; 139:306; 307:384};
period = {[0 99999];[0 99999];[0 99999];[0 99999];[0 99999];[0 99999];[0 99999];[0 99999]};

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


db = commonDB(db, animal, shanks, shank, V1shank, shankGroups, chN, LFPch, eyeCameraCh, chOI, period);


for i = 1:length(db)  
  if ~isfield(db(i), 'series') || isempty(db(i).series)
    db(i).series = 1;
  end
  if ~isfield(db(i), 'evtFilename') || isempty(db(i).evtFilename)
    db(i).evtFilename = [];
  end
  if ~isfield(db(i), 'selectedIntervals') || isempty(db(i).selectedIntervals)
    db(i).selectedIntervals = [];
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
  if isempty(db(i).shankGroups)
    db(i).shankGroups{1} = db(i).shank;
  end
end


clear animal chN eyeCameraCh i LFPch chOI shank shankGroups shanks V1shank seriesFirst seriesLast
fprintf('Created db\n')
