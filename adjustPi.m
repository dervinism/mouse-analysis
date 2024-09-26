function adjHisto = adjustPi(histo)
% A helper function for globalFigs. It adjusts phase histograms so they
% start at -pi/2 rather than -pi.

adjHisto = [histo(6:end,:); histo(1:5,:)];