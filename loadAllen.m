% Allen
type = 'units';
typeRun = 'noRun';
sourceDir = 'S:\cortical_dynamics\User\md406\runNSG';
targetDir = 'S:\cortical_dynamics\User\md406';

allData.M766640955 = loadDataStruct([sourceDir filesep 'runNSG_M766640955\M766640955.mat'], type, typeRun);
allData.M767871931 = loadDataStruct([sourceDir filesep 'runNSG_M767871931\M767871931.mat'], type, typeRun);
allData.M768515987 = loadDataStruct([sourceDir filesep 'runNSG_M768515987\M768515987.mat'], type, typeRun);
allData.M771160300 = loadDataStruct([sourceDir filesep 'runNSG_M771160300\M771160300.mat'], type, typeRun);
allData.M771990200 = loadDataStruct([sourceDir filesep 'runNSG_M771990200\M771990200.mat'], type, typeRun);
allData.M774875821 = loadDataStruct([sourceDir filesep 'runNSG_M774875821\M774875821.mat'], type, typeRun);
allData.M778240327 = loadDataStruct([sourceDir filesep 'runNSG_M778240327\M778240327.mat'], type, typeRun);
allData.M778998620 = loadDataStruct([sourceDir filesep 'runNSG_M778998620\M778998620.mat'], type, typeRun);
allData.M779839471 = loadDataStruct([sourceDir filesep 'runNSG_M779839471\M779839471.mat'], type, typeRun);
allData.M781842082 = loadDataStruct([sourceDir filesep 'runNSG_M781842082\M781842082.mat'], type, typeRun);
allData.M786091066 = loadDataStruct([sourceDir filesep 'runNSG_M786091066\M786091066.mat'], type, typeRun);
allData.M787025148 = loadDataStruct([sourceDir filesep 'runNSG_M787025148\M787025148.mat'], type, typeRun);
allData.M789848216 = loadDataStruct([sourceDir filesep 'runNSG_M789848216\M789848216.mat'], type, typeRun);
allData.M793224716 = loadDataStruct([sourceDir filesep 'runNSG_M793224716\M793224716.mat'], type, typeRun);
allData.M794812542 = loadDataStruct([sourceDir filesep 'runNSG_M794812542\M794812542.mat'], type, typeRun);
allData.M816200189 = loadDataStruct([sourceDir filesep 'runNSG_M816200189\M816200189.mat'], type, typeRun);
allData.M819186360 = loadDataStruct([sourceDir filesep 'runNSG_M819186360\M819186360.mat'], type, typeRun);
allData.M819701982 = loadDataStruct([sourceDir filesep 'runNSG_M819701982\M819701982.mat'], type, typeRun);
allData.M821695405 = loadDataStruct([sourceDir filesep 'runNSG_M821695405\M821695405.mat'], type, typeRun);
allData.M829720705 = loadDataStruct([sourceDir filesep 'runNSG_M829720705\M829720705.mat'], type, typeRun);
allData.M831882777 = loadDataStruct([sourceDir filesep 'runNSG_M831882777\M831882777.mat'], type, typeRun);
allData.M835479236 = loadDataStruct([sourceDir filesep 'runNSG_M835479236\M835479236.mat'], type, typeRun);
allData.M839068429 = loadDataStruct([sourceDir filesep 'runNSG_M839068429\M839068429.mat'], type, typeRun);
allData.M839557629 = loadDataStruct([sourceDir filesep 'runNSG_M839557629\M839557629.mat'], type, typeRun);
allData.M840012044 = loadDataStruct([sourceDir filesep 'runNSG_M840012044\M840012044.mat'], type, typeRun);
allData.M847657808 = loadDataStruct([sourceDir filesep 'runNSG_M847657808\M847657808.mat'], type, typeRun);

% Save the full data set
if strcmpi(typeRun, 'noRun')
  if strcmpi(type, 'all')
    save([targetDir filesep 'allData_allensdk_noRun.mat'], 'allData', '-v7.3');
  elseif strcmpi(type, 'units')
    save([targetDir filesep 'allData_allensdk_noRun.mat'], 'allData', '-v7.3');
  elseif strcmpi(type, 'muas')
    save([targetDir filesep 'allDataMUAs_allensdk_noRun.mat'], 'allData', '-v7.3');
  end
else
  if strcmpi(type, 'all')
    save([targetDir filesep 'allData_allensdk.mat'], 'allData', '-v7.3');
  elseif strcmpi(type, 'units')
    save([targetDir filesep 'allData_allensdk.mat'], 'allData', '-v7.3');
  elseif strcmpi(type, 'muas')
    save([targetDir filesep 'allDataMUAs_allensdk.mat'], 'allData', '-v7.3');
  end
end