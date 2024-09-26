% This script combines two files containing incomplete AnPSD data

fileContents1 = load('filename1');
fileContents2 = load('filename2');

dataStruct.seriesData = fileContents1.dataStruct.seriesData;
dataStruct.popData = fileContents1.dataStruct.popData;
dataStruct.seriesData_ca = fileContents2.dataStruct.popData;

save('filename3','dataStruct','-v7.3');