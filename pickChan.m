function iCh = pickChan(area, animal, chOI)
% iCh = pickChan(area, animal, chOI)
% pickChan function returns LFP channel number given area, animal, and
% corresponding channels of interest.

iCh = [];
if area == 1 || area >= 7
  iCh = round(numel(chOI)/2);
elseif area == 2
  iCh = 1;
elseif area == 4
  if strcmpi(animal, 'M190128_B_MD') || strcmpi(animal, 'M190322_B_MD') || strcmpi(animal, 'M190503_B_MD')
    iCh = 2;
  else
    iCh = 1;
  end
elseif area == 5
  if strcmpi(animal, 'M190523_A_MD')
    iCh = 2;
  else
    iCh = 1;
  end
end