function fileList = addParpoolDependencies(poolObj, path2add, exceptions)

fileList = {};
path2add = regexp(path2add, pathsep, 'split');
not2add = contains(path2add,exceptions);
path2add = path2add(~not2add);
for iPath = 1:numel(path2add)
  folderContents = dir(path2add{iPath});
  for iCont = 1:numel(folderContents)
    fileOrFolder = folderContents(iCont).name;
    if endsWith(fileOrFolder, '.m')
      addAttachedFiles(poolObj,fileOrFolder);
      fileList{numel(fileList)+1} = fileOrFolder; %#ok<*AGROW>
    end
  end
end