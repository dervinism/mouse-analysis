function displayVideoFrames(filename, frameRange)

% Read frames
v = VideoReader(filename);
if numel(frameRange) == 1
  frameRange = [frameRange frameRange];
end
frames = read(v, frameRange);

% Display frames
for frame = 1:numel(frameRange(1):frameRange(2))
  figure;
  if numel(frameRange) == 1
    imagesc(frames)
  else
    imagesc(frames(:,:,:,frame))
  end
end