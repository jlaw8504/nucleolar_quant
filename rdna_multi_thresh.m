function [otsu_area_array, area_mat, thresh_array, otsu_thresh_array] = rdna_multi_thresh(directory)
%% Read in image
gfp_files = dir(fullfile(directory, '*.tif'));
otsu_area_array = zeros([1,size(gfp_files,1)]);
thresh_array = 0:0.01:1;
area_mat = zeros([size(gfp_files,1),length(thresh_array)]);
otsu_thresh_array = zeros([size(gfp_files,1), 1]);
for n = 1:size(gfp_files,1)
    filename = fullfile(gfp_files(n).folder,gfp_files(n).name);
    %open image with readTiffStack
    im = readTiffStack(filename);
    %convert image to double
    imdbl = double(im);
    %% Generate MIP
    %create a max-projection of the image (mip)
    mip = max(imdbl,[],3);
    %% Find and NaN out zero values from MetaMorph Cropping
    if sum(mip(:) == 0) > 0
        mip(mip == 0) = nan;
        warning('%s contains values of zeros', filename);
    end
    %% Normalize image
    %set min pixel to zero
    sub_mip = mip - min(mip(:));
    %set max pixel to one
    norm_mip = sub_mip/max(sub_mip(:));
    %% OTSU thresholding
    otsu_thresh_array(n) = multithresh(norm_mip);
    otsu_bin = norm_mip(norm_mip >= otsu_thresh_array(n));
    otsu_area_array(n) = numel(otsu_bin);
    %% Calc %pixel below threshold at multiple thresholds
    for i = 1:size(thresh_array,2)
        area_mat(n,i) = numel(norm_mip(norm_mip >= thresh_array(i)));
    end
end