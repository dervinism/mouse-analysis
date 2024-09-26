function positions = unitPos(dirname, shCh, probeFile)
% A helper function to loadAsMUA_noResClu. Finds unit channel locations
% including MUAs.

sp = loadKSdir(dirname);

clu = sp.clu; %readNPY([dirname filesep 'spike_clusters.npy']);
res = sp.st * sp.sample_rate; %readNPY([dirname filesep   'spike_times.npy']); 
tmpl = sp.spikeTemplates+1;% readNPY([dirname filesep   'spike_templates.npy']); tmpl = tmpl+1; % they start from 0 (python way)
templateWaveforms = sp.temps; %readNPY([dirname filesep   'templates.npy']); % templates x time x channel
ycoordsCh = sp.ycoords;
xcoordsCh = sp.xcoords;

if isempty(clu) && isempty(res) && isempty(tmpl)
  positions = [];
  return
else
  assert(numel(clu) == numel(res) && numel(res) == numel(tmpl) && max(tmpl) <= size(templateWaveforms, 1))
end

cids = sp.cids; cgs = sp.cgs; %[cids, cgs] = readClusterGroupsCSV([dirname filesep 'cluster_groups.csv']);
%[cids, cgs] = readClusterGroupsCSV([dirname filesep 'cluster_group.tsv']);

uClu = double(unique(clu));

assert(max(abs(double(sort(torow(uClu))) - double(sort(torow(cids))))) == 0, 'should be fully compatible')
assert(~any(cgs >= 3), 'unsorted units remain')

% Make sure no unit is named 0 or 1
if any(uClu == 0)
  m = max(uClu) + 1;
  clu(clu == 0) = m;
  cids(cids == 0) = m;
  uClu(uClu == 0) = m;
end
if any(uClu == 1)
  m = max(uClu) + 1;
  clu(clu == 1) = m;
  cids(cids == 1) = m;
  uClu(uClu == 1) = m;
end
templateWaveforms2D = reshape(templateWaveforms, size(templateWaveforms, 1), []);

load(probeFile, 'ycoords','xcoords')
ycoordsUnique = unique(ycoords);
ycoordsCount = zeros(size(ycoords));
for iCoord = 1:numel(ycoordsUnique)
  ycoordsCount(ycoords == ycoordsUnique(iCoord)) = sum(ycoords == ycoordsUnique(iCoord));
end
if numel(ycoords) <= 65
  ycoordsCount = ycoordsCount./min(ycoordsCount(ycoordsCount > 0));
end

positions = zeros(numel(uClu),4);
uCount = 0;
for u = torow(uClu)
  uCount = uCount + 1;
  h = histc(tmpl(clu == u), 1:size(templateWaveforms, 1)); h = h/sum(h); %#ok<*HISTC>
  w = squeeze(reshape(torow(h)*templateWaveforms2D, size(templateWaveforms, 2), size(templateWaveforms, 3)));
  % approximation for the average waveform of this unit. The computation is
  % done this way because some units can be assigned to more than one
  % template (e.g. after merge in phy)
  
  [~, pos] = max(abs(w(:))); if numel(pos) > 1; pos = pos(1); end
  pos = ceil(pos / size(w, 1));
  ycoordCh = ycoordsCh(pos);
  if exist('ycoordsCount','var')
    ycoordChCount = ycoordsCount(pos);
  else
    ycoordChCount = 1;
  end
  xcoordCh = xcoordsCh(pos);
  posY = find(ycoordCh == ycoords);
  posX = find(xcoordCh == xcoords(posY));
  pos = posY(1) + floor((posX - 1)/ycoordChCount)*shCh;
  positions(uCount,1:3) = [u, pos, pos-1];
  if cgs(cids == u) == 0     % it's noise
    positions(uCount,4) = 0;
  elseif cgs(cids == u) == 1 % it's MUA
    positions(uCount,4) = 1;
  else
    positions(uCount,4) = u;
  end
end