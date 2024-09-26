function shiftedValue = shiftRadAxis(originalValue, axisLims)
% A function for rotating radian axis values. It is a helper function of
% halfCorrPlot.

assert(abs(axisLims(2)-axisLims(1)) == 2*pi, 'Axis limits are not exactly 2 pi');

if originalValue < axisLims(1)
  shiftedValue = axisLims(2) - (axisLims(1) - originalValue);
elseif originalValue > axisLims(2)
  shiftedValue = axisLims(1) + (originalValue - axisLims(2));
else
  shiftedValue = originalValue;
end