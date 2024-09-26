poolObj = parpool(24);
%poolObj = gcp;

if exist('path2add','var')
  exceptions = {'chronux','matlib'};
  filesAttached = addParpoolDependencies(poolObj, path2add, exceptions);
end