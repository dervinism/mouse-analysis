% This function removes excess lfp channel analysis data caused by a bug in
% the lfpLoad function.

dataFile = 'M190114_A_MD';

load(dataFile);

fnsData = fieldnames(dataStruct.seriesData);
for dbCount = 1:numel(fnsData)
  dbStruct = dataStruct.seriesData.(fnsData{dbCount});
  entryName = dbStruct.db(dbCount).entryName;
  
  if isfield(dbStruct, 'lfpPowerData')
    chOIDB = dbStruct.lfpPowerData.chOIDB;
    rippleRateTemp = dbStruct.lfpPowerData.rippleRate;
    meanRippleRateTemp = dbStruct.lfpPowerData.meanRippleRate;
    theta2deltaRatioTemp = dbStruct.lfpPowerData.theta2deltaRatio;
    slowPowerTemp = dbStruct.lfpPowerData.slowPower;
    fastPowerTemp = dbStruct.lfpPowerData.fastPower;
    ultraFastPowerTemp = dbStruct.lfpPowerData.ultraFastPower;
    dbStruct.lfpPowerData.rippleRate = {};
    dbStruct.lfpPowerData.meanRippleRate = {};
    dbStruct.lfpPowerData.theta2deltaRatio = {};
    dbStruct.lfpPowerData.slowPower = {};
    dbStruct.lfpPowerData.fastPower = {};
    dbStruct.lfpPowerData.ultraFastPower = {};
    for iCh = 1:numel(chOIDB)
      dbStruct.lfpPowerData.rippleRate{iCh} = rippleRateTemp{iCh};
      dbStruct.lfpPowerData.meanRippleRate{iCh} =  meanRippleRateTemp{iCh};
      dbStruct.lfpPowerData.theta2deltaRatio{iCh} = theta2deltaRatioTemp{iCh};
      dbStruct.lfpPowerData.slowPower{iCh} = slowPowerTemp{iCh};
      dbStruct.lfpPowerData.fastPower{iCh} = fastPowerTemp{iCh};
      dbStruct.lfpPowerData.ultraFastPower{iCh} = ultraFastPowerTemp{iCh};
    end
  else
    continue
  end
  
  if isfield(dbStruct, 'lfpphaseCohDataMotion')
    rippleRateTemp = dbStruct.lfpphaseCohDataMotion.rippleRate;
    t2dRatioTemp = dbStruct.lfpphaseCohDataMotion.t2dRatio;
    slowPowerTemp = dbStruct.lfpphaseCohDataMotion.slowPower;
    fastPowerTemp = dbStruct.lfpphaseCohDataMotion.fastPower;
    ultraFastPowerTemp = dbStruct.lfpphaseCohDataMotion.ultraFastPower;
    dbStruct.lfpphaseCohDataMotion.rippleRate = {};
    dbStruct.lfpphaseCohDataMotion.t2dRatio = {};
    dbStruct.lfpphaseCohDataMotion.slowPower = {};
    dbStruct.lfpphaseCohDataMotion.fastPower = {};
    dbStruct.lfpphaseCohDataMotion.ultraFastPower = {};
    for iCh = 1:numel(chOIDB)
      dbStruct.lfpphaseCohDataMotion.rippleRate{iCh} = rippleRateTemp{iCh};
      dbStruct.lfpphaseCohDataMotion.t2dRatio{iCh} = t2dRatioTemp{iCh};
      dbStruct.lfpphaseCohDataMotion.slowPower{iCh} = slowPowerTemp{iCh};
      dbStruct.lfpphaseCohDataMotion.fastPower{iCh} = fastPowerTemp{iCh};
      dbStruct.lfpphaseCohDataMotion.ultraFastPower{iCh} = ultraFastPowerTemp{iCh};
    end
  end
  
  if isfield(dbStruct, 'lfpphaseCohData')
    rippleRateTemp = dbStruct.lfpphaseCohData.rippleRate;
    t2dRatioTemp = dbStruct.lfpphaseCohData.t2dRatio;
    slowPowerTemp = dbStruct.lfpphaseCohData.slowPower;
    fastPowerTemp = dbStruct.lfpphaseCohData.fastPower;
    ultraFastPowerTemp = dbStruct.lfpphaseCohData.ultraFastPower;
    dbStruct.lfpphaseCohData.rippleRate = {};
    dbStruct.lfpphaseCohData.t2dRatio = {};
    dbStruct.lfpphaseCohData.slowPower = {};
    dbStruct.lfpphaseCohData.fastPower = {};
    dbStruct.lfpphaseCohData.ultraFastPower = {};
    for iCh = 1:numel(chOIDB)
      dbStruct.lfpphaseCohData.rippleRate{iCh} = rippleRateTemp{iCh};
      dbStruct.lfpphaseCohData.t2dRatio{iCh} = t2dRatioTemp{iCh};
      dbStruct.lfpphaseCohData.slowPower{iCh} = slowPowerTemp{iCh};
      dbStruct.lfpphaseCohData.fastPower{iCh} = fastPowerTemp{iCh};
      dbStruct.lfpphaseCohData.ultraFastPower{iCh} = ultraFastPowerTemp{iCh};
    end
  end
  
  if isfield(dbStruct, 'lfpphaseCohDataPR')
    rippleRateTemp = dbStruct.lfpphaseCohDataPR.rippleRate;
    t2dRatioTemp = dbStruct.lfpphaseCohDataPR.t2dRatio;
    slowPowerTemp = dbStruct.lfpphaseCohDataPR.slowPower;
    fastPowerTemp = dbStruct.lfpphaseCohDataPR.fastPower;
    ultraFastPowerTemp = dbStruct.lfpphaseCohDataPR.ultraFastPower;
    dbStruct.lfpphaseCohDataPR.rippleRate = {};
    dbStruct.lfpphaseCohDataPR.t2dRatio = {};
    dbStruct.lfpphaseCohDataPR.slowPower = {};
    dbStruct.lfpphaseCohDataPR.fastPower = {};
    dbStruct.lfpphaseCohDataPR.ultraFastPower = {};
    for iCh = 1:numel(chOIDB)
      dbStruct.lfpphaseCohDataPR.rippleRate{iCh} = rippleRateTemp{iCh};
      dbStruct.lfpphaseCohDataPR.t2dRatio{iCh} = t2dRatioTemp{iCh};
      dbStruct.lfpphaseCohDataPR.slowPower{iCh} = slowPowerTemp{iCh};
      dbStruct.lfpphaseCohDataPR.fastPower{iCh} = fastPowerTemp{iCh};
      dbStruct.lfpphaseCohDataPR.ultraFastPower{iCh} = ultraFastPowerTemp{iCh};
    end
  end
  
  dataStruct.seriesData.(fnsData{dbCount}) = dbStruct;
end

save(dataFile,'dataStruct','-v7.3');