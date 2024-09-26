function db = commonDB(db, animal, repository, researcher, shanks, shank, chN, LFPch, eyeCameraCh, chOI, period, periodLFP)
% commonDB is a helper function for makedb. It assigns db fields that are
% shared across different recording sessions in the same animal.


for i = 1:numel(db)
  db(i).animal = animal;
  db(i).repository = repository;
  db(i).researcher = researcher;
  db(i).shanks = shanks(i);  % Number of shanks
  db(i).shank = shank{i};    % Shank IDs
  db(i).chN = chN{i};
  db(i).LFPch = LFPch{i};
  db(i).eyeCameraCh = eyeCameraCh{i};
  db(i).chOI = chOI{i};
  db(i).period = period{i};
  db(i).periodLFP = periodLFP{i};
end