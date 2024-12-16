type = 'units';
sourceDirNeuronexus = 'D:\infraslow-dynamics\04_data_analysis\001_uol';
sourceDirNeuropixels = 'D:\infraslow-dynamics\04_data_analysis\001_uol';
targetDir = 'D:\infraslow-dynamics\04_data_analysis\001_uol\noRun';

% UOL Neuronexus
allData.M180313 = loadDataStruct([sourceDirNeuronexus filesep 'M180313\M180313.mat'], type);
allData.M180711_MD = loadDataStruct([sourceDirNeuronexus filesep 'M180711_MD\M180711_MD.mat'], type);
allData.M180712_MD = loadDataStruct([sourceDirNeuronexus filesep 'M180712_MD\M180712_MD.mat'], type);
allData.M180713_MD = loadDataStruct([sourceDirNeuronexus filesep 'M180713_MD\M180713_MD.mat'], type);
allData.M180717_B_MD = loadDataStruct([sourceDirNeuronexus filesep 'M180717_B_MD\M180717_B_MD.mat'], type);
allData.M180719_MD = loadDataStruct([sourceDirNeuronexus filesep 'M180719_MD\M180719_MD.mat'], type);
allData.M180924_A_MD = loadDataStruct([sourceDirNeuronexus filesep 'M180924_A_MD\M180924_A_MD.mat'], type);
allData.M180924_B_MD = loadDataStruct([sourceDirNeuronexus filesep 'M180924_B_MD\M180924_B_MD.mat'], type);
allData.M181112_B_MD = loadDataStruct([sourceDirNeuronexus filesep 'M181112_B_MD\M181112_B_MD.mat'], type);
allData.M181210_B_MD = loadDataStruct([sourceDirNeuronexus filesep 'M181210_B_MD\M181210_B_MD.mat'], type);
allData.M190114_A_MD = loadDataStruct([sourceDirNeuronexus filesep 'M190114_A_MD\M190114_A_MD.mat'], type);
allData.M190128_A_MD = loadDataStruct([sourceDirNeuronexus filesep 'M190128_A_MD\M190128_A_MD.mat'], type);
allData.M190128_B_MD = loadDataStruct([sourceDirNeuronexus filesep 'M190128_B_MD\M190128_B_MD.mat'], type);
allData.M190128_C_MD = loadDataStruct([sourceDirNeuronexus filesep 'M190128_C_MD\M190128_C_MD.mat'], type);
allData.M190322_B_MD = loadDataStruct([sourceDirNeuronexus filesep 'M190322_B_MD\M190322_B_MD.mat'], type);
allData.M190322_C_MD = loadDataStruct([sourceDirNeuronexus filesep 'M190322_C_MD\M190322_C_MD.mat'], type);
allData.M190503_A_MD = loadDataStruct([sourceDirNeuronexus filesep 'M190503_A_MD\M190503_A_MD.mat'], type);
allData.M190503_B_MD = loadDataStruct([sourceDirNeuronexus filesep 'M190503_B_MD\M190503_B_MD.mat'], type);
allData.M190523_A_MD = loadDataStruct([sourceDirNeuronexus filesep 'M190523_A_MD\M190523_A_MD.mat'], type);
allData.M190523_B_MD = loadDataStruct([sourceDirNeuronexus filesep 'M190523_B_MD\M190523_B_MD.mat'], type);

% UOL Neuropixels
allData.M191018_MD = loadDataStruct([sourceDirNeuropixels filesep 'M191018_MD\M191018_MD.mat'], type);
allData.M191106_MD = loadDataStruct([sourceDirNeuropixels filesep 'M191106_MD\M191106_MD.mat'], type);
allData.M191107_MD = loadDataStruct([sourceDirNeuropixels filesep 'M191107_MD\M191107_MD.mat'], type);
allData.M191119_A_MD = loadDataStruct([sourceDirNeuropixels filesep 'M191119_A_MD\M191119_A_MD.mat'], type);
allData.M191119_B_MD = loadDataStruct([sourceDirNeuropixels filesep 'M191119_B_MD\M191119_B_MD.mat'], type);
allData.M191119_C_MD = loadDataStruct([sourceDirNeuropixels filesep 'M191119_C_MD\M191119_C_MD.mat'], type);
allData.M191128_A_MD = loadDataStruct([sourceDirNeuropixels filesep 'M191128_A_MD\M191128_A_MD.mat'], type);
allData.M191128_B_MD = loadDataStruct([sourceDirNeuropixels filesep 'M191128_B_MD\M191128_B_MD.mat'], type);
allData.M200316_MD = loadDataStruct([sourceDirNeuropixels filesep 'M200316_MD\M200316_MD.mat'], type);
allData.M200317_MD = loadDataStruct([sourceDirNeuropixels filesep 'M200317_MD\M200317_MD.mat'], type);
allData.M200318_MD = loadDataStruct([sourceDirNeuropixels filesep 'M200318_MD\M200318_MD.mat'], type);
allData.M200319_MD = loadDataStruct([sourceDirNeuropixels filesep 'M200319_MD\M200319_MD.mat'], type);
allData.M200323_MD = loadDataStruct([sourceDirNeuropixels filesep 'M200323_MD\M200323_MD.mat'], type);
allData.M200324_MD = loadDataStruct([sourceDirNeuropixels filesep 'M200324_MD\M200324_MD.mat'], type);
allData.M200325_MD = loadDataStruct([sourceDirNeuropixels filesep 'M200325_MD\M200325_MD.mat'], type);
allData.M210802_MD = loadDataStruct([sourceDirNeuropixels filesep 'M210802_MD\M210802_MD.mat'], type);

% Save the full data set
if strcmpi(type, 'all')
  save([targetDir filesep 'allData_uol.mat'], 'allData', '-v7.3'); %#ok<*UNRCH>
elseif strcmpi(type, 'units')
  save([targetDir filesep 'allData_uol.mat'], 'allData', '-v7.3');
elseif strcmpi(type, 'muas')
  save([targetDir filesep 'allDataMUAs_uol.mat'], 'allData', '-v7.3');
end