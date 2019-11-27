function corr_array = rdna_corr(directory, rdna_pattern, nuc_pattern)
%%corr_array Calculates Pearson's linear correlation coefficient for
%%fluorescent signals above a background.
%
%   input :
%       directory : A string variable pointing to the directory containing
%       the cropped, fluorescent images.
%
%       rdna_pattern : A string variable that contains a pattern that points
%       to a particular image channel that contains the rDNA fluorescent
%       signal (from the lacO/LacI-GFP rDNA array). For example 'GFP.tif'.
%
%       nuc_pattern : A string variable that contains a pattern that points
%       to a particular image channel containing the signal of the
%       fluorescent protein you wish to compare to the rDNA signal.
%
%   output :
%       corr_array : An array variable containing the correlation
%       coefficients of each image comparison.
%% Parse File directory
rdna_files = dir(fullfile(directory, rdna_pattern));
nuc_files = dir(fullfile(directory, nuc_pattern));
%% Size Check
if numel(rdna_files) ~= numel(nuc_files)
    error('Different number of images for rDNA locus and nucleolar protein');
end
%% Pre-allocate corr array
corr_array = zeros([1, numel(rdna_files)]);
%% Loop over the files of Max Intensity Projections
for n = 1:numel(rdna_files)
    rdna_mat = readTiffStack(fullfile(directory, rdna_files(n).name));
    nuc_mat = readTiffStack(fullfile(directory, nuc_files(n).name));
    %% Create binary index from rDNA image
    rdna_bin = rdna_mat > multithresh(rdna_mat);
    %% Use rdna_bin to create list of intensity values
    rdna_overlap = rdna_mat(rdna_bin);
    nuc_overlap = nuc_mat(rdna_bin);
    %% Use corr to calculate Pearson's linear correlation coefficient
    rho = corr(rdna_overlap, nuc_overlap);
    corr_array(n) = rho;
end