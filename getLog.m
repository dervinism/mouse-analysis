function logData = getLog(data)

logData = log10(data);
logData(isinf(logData)) = -3;
end