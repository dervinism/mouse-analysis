function sz = h5_dataset_size(h5_path, dataset_path)
% H5_DATASET_SIZE  Get a dataset's dimensions via low-level HDF5 (fast).
%
%   sz = h5_dataset_size(h5_path, dataset_path)
%
%   Low-level replacement for h5info(h5_path, dataset_path).Dataspace.Size.
%   MATLAB's high-level h5info appears to eagerly inspect a dataset's full
%   storage/chunk layout, which can take minutes for huge, finely chunked
%   arrays nested inside large v7.3 MAT-files. This uses H5D.open +
%   H5D.get_space + H5S.get_simple_extent_dims, which only reads the
%   dataspace descriptor.
%
%   Dimensions are returned in the same order h5info reports them (MATLAB/
%   column-major convention), i.e. reversed from HDF5's own C-order storage
%   dims -- this is a drop-in replacement for existing
%   h5info(...).Dataspace.Size call sites and matches size(...) semantics
%   (sz(2) is the number of columns of a 2-D dataset, etc.).
%
%   Inputs:
%       h5_path      - Path to the HDF5 file (a v7.3 MAT-file is a valid
%                      HDF5 file and can be used directly here).
%       dataset_path - Absolute path of the dataset, e.g.
%                      '/allData/M180313/dataStruct/seriesData/<series>/shankData/shank1/MUAs'.
%
%   Output:
%       sz - Row vector of dimension sizes, h5info-Size-/size()-compatible order.
%
%   Example:
%       sz = h5_dataset_size(matFile, [seriesPath '/shankData/shank1/MUAs']);
%       nSamples = sz(2);

    file_id = H5F.open(h5_path, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
    cleanupFile = onCleanup(@() H5F.close(file_id)); %#ok<NASGU>

    dset_id = H5D.open(file_id, dataset_path);
    cleanupDset = onCleanup(@() H5D.close(dset_id)); %#ok<NASGU>

    space_id = H5D.get_space(dset_id);
    cleanupSpace = onCleanup(@() H5S.close(space_id)); %#ok<NASGU>

    [~, dims, ~] = H5S.get_simple_extent_dims(space_id);
    sz = fliplr(dims);
end
