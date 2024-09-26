% Computes total movement (using simple derivative and some fancy MATLAB algo)

clear all %#ok<CLALL>

eyeData = {'R:\Neuropix\Shared\Data\M200325_MD\20200325122741\video_52397_bad_quality.mj2'};
         
for dbCount = 1:numel(eyeData)
  v = VideoReader(eyeData{dbCount}); %#ok<TNMLP>
  [filepath,name,ext] = fileparts(eyeData{dbCount});
  resultsFilename = [filepath filesep name '_tm.mat'];
  
  % This is Michael's code onwards
  if exist(resultsFilename, 'file')    
    load(resultsFilename);
  end
  if exist('s', 'var') && exist('sa', 'var')  
    fprintf('- already done\n')
    clearvars -except eyeData dbCount
    continue % stuff for this video already computed
  end
  if ~exist('s', 'var')
    foregroundDetector = vision.ForegroundDetector();
    s = zeros(1, v.NumberOfFrames); %#ok<*VIDREAD>
  else
    foregroundDetector = [];
  end
  sa = zeros(1, v.NumberOfFrames);
  y = v.read(1);
  for i = 1:v.NumberOfFrames
    x = v.read(i);  
    sa(i) = sum(abs(x(:)-y(:))); % total difference wrt previous frame
    y = x;
    if ~isempty(foregroundDetector)
      fg = foregroundDetector.step(x);
      s(i) = sum(fg(:));
    end
    if mod(i, 5000) == 0
      disp([num2str(i) ' / ' num2str(v.NumberOfFrames)]);
    end
  end
  
  % Z-score
  s = (s - mean(s))/std(s); % total movement by using vision.ForegroundDetector
  sa = (sa - mean(sa))/std(sa); % total movement by simple derivative (diff)
  
  clear v results fg foregroundDetector
  save(resultsFilename, 's', 'sa'); % "total motion" data
  fprintf('\n')
  clearvars -except eyeData dbCount
end