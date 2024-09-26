 function corrChannels(filename, nChans, chans2ignore)
% The function reads binary files and correlates channels.
%
% Input: filename should include the full path and the extension.
%        nChans is the total number of channels in the recording.
%        chans2ignore tells which channels should be excluded from analysis.
%          Typically these are extra channels that do not correspond to
%          probe LFP recording channels.
%
% Output: figs contains figure handles.

chunkSize = 10000000;

fid = [];

d = dir(filename);
nSampsTotal = d.bytes/nChans/2;
nChunksTotal = ceil(nSampsTotal/chunkSize);

chans = ones(size(1:nChans));
chans(chans2ignore) = 0;
chans = logical(chans);
pairs = [(1:sum(chans))' ([2:sum(chans) 1])'];
corrs = [pairs zeros(size(pairs,1),1)];

filepath = fileparts(filename);

try
  
  fid = fopen(filename, 'r');
  
  chunkInd = 1;
  while 1
    fprintf(1, 'chunk %d/%d\n', chunkInd, nChunksTotal);
    dat = fread(fid, [nChans chunkSize], '*int16');
    if ~isempty(dat)
      dat = dat(chans,:);
      
      % Correlate channels
      pairs = [(1:sum(chans))' ([2:sum(chans) 1])'];
      corrCoefs = zeros(size(pairs,1),1);
      for p = 1:size(pairs,1)
        fprintf(1, 'pair %d/%d\n', p, size(pairs,1));
        corrCoefs(p) = corrSimple(double(dat(pairs(p,1),:)), double(dat(pairs(p,2),:)), 'Pearson');
      end
      corrs(:,3) = corrs(:,3) + corrCoefs;
      
    else
      break
    end
    chunkInd = chunkInd+1;
  end
  corrs(:,3) = corrs(:,3)./(chunkInd-1); %#ok<NASGU>
  
  save([filepath filesep 'chanCorrs.mat'], 'corrs');
  
  close(gcf)
  fclose(fid);
  
catch me
  
  if ~isempty(fid)
    fclose(fid);
  end
  
  rethrow(me)
  
end
