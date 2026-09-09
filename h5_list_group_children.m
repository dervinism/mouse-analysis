function names = h5_list_group_children(h5_path, group_path)
% H5_LIST_GROUP_CHILDREN  List immediate child names of an HDF5 group (fast).
%
%   names = h5_list_group_children(h5_path, group_path)
%
%   Low-level replacement for {h5info(h5_path, group_path).Groups.Name} or
%   .Datasets.Name. MATLAB's high-level h5info appears to eagerly inspect the
%   full storage/chunk layout of every descendant object, which can take
%   minutes when a group contains huge datasets nested many levels deep (as
%   v7.3 MAT-files holding large struct trees do) -- even when only the
%   immediate child names are needed. This uses H5G.get_info +
%   H5L.get_name_by_idx, which only reads the group's link table and is
%   near-instant regardless of how large any descendant dataset is.
%
%   Inputs:
%       h5_path    - Path to the HDF5 file (a v7.3 MAT-file is a valid
%                    HDF5 file and can be used directly here).
%       group_path - Absolute path of the group to list, e.g. '/allData'
%                    or '/allData/M180313/dataStruct/seriesData'.
%
%   Output:
%       names - Cell array of bare child names (e.g. 'M180313', not
%               '/allData/M180313'), in link-creation order, regardless of
%               whether each child is itself a group or a dataset. Prefix
%               with group_path yourself if a full path is needed (e.g.
%               [group_path '/' names{i}]).
%
%   Example:
%       animalIDs = h5_list_group_children(matFile, '/allData');
%       seriesNames = h5_list_group_children(matFile, ...
%           ['/allData/' animalIDs{1} '/dataStruct/seriesData']);

    file_id = H5F.open(h5_path, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
    cleanupFile = onCleanup(@() H5F.close(file_id)); %#ok<NASGU>

    grp_id = H5G.open(file_id, group_path);
    cleanupGrp = onCleanup(@() H5G.close(grp_id)); %#ok<NASGU>

    ginfo = H5G.get_info(grp_id);
    names = cell(ginfo.nlinks, 1);
    for i = 0:ginfo.nlinks-1
        names{i+1} = H5L.get_name_by_idx(file_id, group_path, ...
            'H5_INDEX_NAME', 'H5_ITER_INC', i, 'H5P_DEFAULT');
    end
end
