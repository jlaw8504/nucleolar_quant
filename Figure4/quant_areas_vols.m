function [infocus_areas, mip_areas, sip_areas, vols, int_ints, mean_ints] = quant_areas_vols(pattern)
%%quant_areas_vols calculate the areas in pixels and volume in voxels of 
%%nucleolar signals from cropped cell image stacks
%
%   inputs :
%       pattern : A string variable that specfies which file in the current
%       working directory to parse.
%
%   outputs :
%       infocus_areas = A row vector containing the areas of the signals of
%       the in-focus planes
%
%       mip_areas : A row vector containing the areas of the signals of the
%       maximum intensity projection
%
%       sip_areas : A row vector containing the areas of the signals of the
%       sum intensity projection
%
%       vols : A row vector containing the volume (in voxels) of the entire
%       image stack
%
%       int_ints : A row vector containing the integrated intensity of the
%       nucleolar region in AU
%
%       mean_ints : A row vector containing the mean intensity of all the
%       nucleolar signals after thresholding
%
%   NOTE: Image stacks that have the brightest pixel in either the first or
%   last Z-plane area excluded and their values are set to nan in the
%   output arrays.

rfp_files = dir(pattern);
infocus_areas = zeros([numel(rfp_files), 1]);
mip_areas = infocus_areas;
sip_areas = infocus_areas;
vols = infocus_areas;
int_ints = infocus_areas;
mean_ints = infocus_areas;

for n = 1:numel(rfp_files)
im_mat = readTiffStack(rfp_files(n).name);
[~, idx] = max(im_mat(:));
[~, ~, z] = ind2sub(size(im_mat), idx);
if z == 1 || z == size(im_mat,3)
    infocus_areas(n) = nan;
    mip_areas(n) = nan;
    sip_areas(n) = nan;
    vols(n) = nan;
    int_ints(n) = nan;
    mean_ints(n) = nan;
    continue
end
infocus = im_mat(:, :, z);
mip = max(im_mat, [], 3);
sip = sum(im_mat, 3);

thresh_array_infocus = multithresh(infocus, 2);
thresh_array_mip = multithresh(mip, 2);
thresh_array_sip = multithresh(sip, 2);

infocus_areas(n) = sum(infocus(:) > thresh_array_infocus(end));
sip_areas(n) = sum(sip(:) > thresh_array_sip(end));
mip_areas(n) = sum(mip(:) > thresh_array_mip(end));

thresh_array = multithresh(im_mat, 2);
bin_mat = im_mat > thresh_array(end);
vols(n) = sum(bin_mat(:));
int_ints(n) = sum(im_mat(bin_mat));
mean_ints(n) = int_ints(n)/vols(n) - thresh_array(1);
end