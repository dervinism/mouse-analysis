function colourCode = areaColours(area)
% colourCode = animalColours(animal) assigns a colour code to an area.
% Input: area - an area number or string name.
% Output: colourCode - a colour code vector.

if (isnumeric(area) && area == 1) || strcmpi(area, 'S1')
  colourCode = [0, 0, 1];
elseif (isnumeric(area) && area == 2) || (~isnumeric(area) && strcmpi(area, 'VB'))
  colourCode = [0, 0.5, 0];
elseif (isnumeric(area) && area == 3) || (~isnumeric(area) && strcmpi(area, 'Th'))
  colourCode = [1, 0, 0];
elseif (isnumeric(area) && area == 4) || (~isnumeric(area) && strcmpi(area, 'DG'))
  colourCode = [0, 0.75, 0.75];
elseif (isnumeric(area) && area == 5) || (~isnumeric(area) && strcmpi(area, 'CA1'))
  colourCode = [0.75, 0, 0.75];
elseif (isnumeric(area) && area == 6) || (~isnumeric(area) && strcmpi(area, 'Hp'))
  colourCode = [0.75, 0.75, 0];
elseif (isnumeric(area) && area == 7) || (~isnumeric(area) && strcmpi(area, 'RSC'))
  colourCode = [0.25, 0.25, 0.25];
elseif (isnumeric(area) && area == 8) || (~isnumeric(area) && strcmpi(area, 'mPFC'))
  colourCode = [0, 0.447, 0.741];
elseif (isnumeric(area) && area == 9) || (~isnumeric(area) && strcmpi(area, 'V1'))
  colourCode = [0.85, 0.325, 0.098];
elseif (isnumeric(area) && area == 10) || (~isnumeric(area) && strcmpi(area, 'Cx'))
  colourCode = [0.929, 0.694, 0.125];
elseif (isnumeric(area) && area == 11) || (~isnumeric(area) && strcmpi(area, 'reserved'))
  colourCode = [0.494, 0.184, 0.556];
elseif (isnumeric(area) && area == 12) || (~isnumeric(area) && strcmpi(area, 'reserved'))
  colourCode = [0.466, 0.674, 0.188];
elseif (isnumeric(area) && area == 13) || (~isnumeric(area) && strcmpi(area, 'reserved'))
  colourCode = [0.301, 0.745, 0.933];
elseif (isnumeric(area) && area == 14) || (~isnumeric(area) && strcmpi(area, 'reserved'))
  colourCode = [0.635, 0.078, 0.184];
elseif (isnumeric(area) && area == 15) || (~isnumeric(area) && strcmpi(area, 'reserved'))
  colourCode = ([0, 0, 1] + [0, 0.5, 0])./2;
elseif (isnumeric(area) && area == 16) || (~isnumeric(area) && strcmpi(area, 'reserved'))
  colourCode = ([1, 0, 0] + [0, 0.75, 0.75])./2;
elseif (isnumeric(area) && area == 17) || (~isnumeric(area) && strcmpi(area, 'reserved'))
  colourCode = ([0.75, 0, 0.75] + [0.75, 0.75, 0])./2;
elseif (isnumeric(area) && area == 18) || (~isnumeric(area) && strcmpi(area, 'reserved'))
  colourCode = ([0.25, 0.25, 0.25] + [0, 0.447, 0.741])./2;
elseif (isnumeric(area) && area == 19) || (~isnumeric(area) && strcmpi(area, 'reserved'))
  colourCode = ([0.85, 0.325, 0.098] + [0.929, 0.694, 0.125])./2;
end