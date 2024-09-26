function [h, phase, envelope] = hilbertTransform(data)
% A helper function of AnHT for obtaining Hilbert transform of series data
% and estimating envelope and phase.

h = hilbert(data');
h = h';
phase = angle(h);
envelope = abs(h);
end