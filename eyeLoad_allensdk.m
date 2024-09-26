function [pupilArea, time] = eyeLoad_allensdk(filename)


%% Load the file with gaze data
load(filename); %#ok<*LOAD>


%% Determine relevant data columns
timeColumn = find(contains(cellstr(headerGaze),'Time (s)'));
pupilAreaColumn = find(contains(cellstr(headerGaze),'raw_pupil_area'));


%% Extract relevant data
time = gazeData(:,timeColumn)'; %#ok<*FNDSB,*USENS>
pupilArea = gazeData(:,pupilAreaColumn)'; % cm^2