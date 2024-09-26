% This function updates data structure with certain recording data
% properties

fclose all;
close all
clear
clc

dataFile = 'M180719_MD';
load(dataFile);
lists;
lfpParams;

fnsData = fieldnames(dataStruct.seriesData);
for dbCount = 1:numel(fnsData)
  if dbCount == 1
    dbStruct = dataStruct.seriesData.(fnsData{dbCount});
    db = dbStruct.db;
    
    for dbCount2 = 1:numel(fnsData)
      entryName = db(dbCount2).entryName;
      seriesName = seriesFromEntry(entryName);
      basefilename = db(dbCount2).basefilename;
      
      % Probe ID and area
      [~, areaID] = determineAreaFromSeries(seriesName);
      [probeID, area] = num2area(areaID);
      db(dbCount2).probeID = probeID;
      db(dbCount2).area = area;
      
      % Condition
      for entry = 1:numel(awake)
        if strcmpi(awake{entry}, seriesName(1:min([14, numel(seriesName)])))
          condition = 'awake';
        end
      end
      for entry = 1:numel(anaesthesia)
        if strcmpi(anaesthesia{entry}, seriesName(1:min([14, numel(seriesName)])))
          condition = 'anaesthesia';
        end
      end
      db(dbCount2).condition = condition;
      
      % LFPOI
      db(dbCount2).LFPOI = chOI{dbCount2};
      
      % Sharp-wave ripples
      for dbCount3 = 1:numel(fullSeries)
        seriesName2 = seriesFromEntry(fullSeries{dbCount3});
        if strcmpi(seriesName, seriesName2)
          db(dbCount2).ripples = 1;
          db(dbCount2).chRipples = chOI{dbCount2}(ch{dbCount3});
          break
        else
          db(dbCount2).ripples = 0;
          db(dbCount2).chRipples = [];
        end
      end
      
      % CAR
      if contains(basefilename, 'NoCAR')
        db(dbCount2).CAR = 0;
      elseif contains(basefilename, 'CAR')
        db(dbCount2).CAR = 1;
      else
        db(dbCount2).CAR = 0;
      end
    end
  end
  
  dataStruct.seriesData.(fnsData{dbCount}).dbSeries = db(dbCount);
  dataStruct.seriesData.(fnsData{dbCount}).db = db;
end

if isfield(dataStruct, 'seriesData_ca')
  fnsData = fieldnames(dataStruct.seriesData_ca);
  for dbCount = 1:numel(fnsData)
    [seriesName1, seriesName2] = seriesNames(fnsData{dbCount});
    [~, areaID1] = determineAreaFromSeries(seriesName1);
    [~, areaID2] = determineAreaFromSeries(seriesName2);
    [probeID1, area1] = num2area(areaID1);
    [probeID2, area2] = num2area(areaID2);
    probe = [num2str(probeID1) 'Vs' num2str(probeID2)];
    area = [area1 'Vs' area2];
    dataStruct.seriesData_ca.(fnsData{dbCount}).probe = probe;
    dataStruct.seriesData_ca.(fnsData{dbCount}).area = area;
  end
end

save(dataFile,'dataStruct','-v7.3');



function [probeID, area] = num2area(areaID)

if strcmpi(areaID, '1')       % S1
  probeID = 1;
  area = 'S1';
elseif strcmpi(areaID, '2')   % VB1
  probeID = 2;
  area = 'VB';
elseif strcmpi(areaID, '3')   % Po
  probeID = 2;
  area = 'Po';
elseif strcmpi(areaID, '4')   % LP1
  probeID = 2;
  area = 'LP';
elseif strcmpi(areaID, '5')   % DG
  probeID = 2;
  area = 'DG';
elseif strcmpi(areaID, '6')   % CA1
  probeID = 2;
  area = 'CA1';
elseif strcmpi(areaID, '7')   % RSC
  probeID = 2;
  area = 'RSC';
elseif strcmpi(areaID, '8')   % VB2
  probeID = 1;
  area = 'VB';
elseif strcmpi(areaID, '9')   % LP2
  probeID = 1;
  area = 'LP';
elseif strcmpi(areaID, '10')  % DLGN
  probeID = 1;
  area = 'DLGN';
elseif strcmpi(areaID, '11')  % CA3
  probeID = 1;
  area = 'CA3';
elseif strcmpi(areaID, '12')  % mPFC
  probeID = 1;
  area = 'mPFC';
elseif strcmpi(areaID, '13')  % lV1
  probeID = 1;
  area = 'lV1';
elseif strcmpi(areaID, '14')  % rV1
  probeID = 2;
  area = 'rV1';
elseif strcmpi(areaID, '24')  % Th1
  probeID = 2;
  area = 'Th';
elseif strcmpi(areaID, '56')  % Hp
  probeID = 2;
  area = 'Hp';
elseif strcmpi(areaID, '810') % Th2
  probeID = 1;
  area = 'Th';
else
  error('Unknown area');
end
end