function [rho, rdna_vol, rdna_std, nuc_vol, nuc_std, dispersion_ratio, nuc_sph] = rnase_quant(rdna_filename, nuc_filename, pixel_size, z_step)
%%noy_quant calculates the volume and standard deviation of image
%%intensities for cells containing both rdna and nucleolar protein signals
%
%   input :
%       rdna_filename : A string variable specifying which rdna image file
%       to parse.
%
%       nuc_filename : A string variable specifying which nucleolar image
%       file to parse.
%
%       pixel_size : A scalar variable specifying the pixel size.
%
%       z_step : A scalar varible specifying the z-step size.
%
%   output :
%       rho : A scalar variable containing Pearson's correlation between
%       the rdna protein signal and the nucleolar protein signal.
%
%       rdna_vol : A scalar variable containing the volume of the rdna
%       protein signal.
%
%       rdna_std : A scalar variable containing the standard deviation of
%       the voxel intensities of the rdna protein signal.
%
%       nuc_vol : A scalar variable containing the volume of the nucleolar
%       protein signal.
%
%       nuc_std : A scalar variable containing the standard deviation of
%       the voxel intensities of the nucleolar protein signal.
%
%       std_displacement : A scalar variable containing the STD of the distance
%       of each rDNA voxel of the largest rDNA signal to the center of the
%       largest nucleolar signal.

%% Parse the Image files
rdna_mat = readTiffStack(rdna_filename);
nuc_mat = readTiffStack(nuc_filename);

%% Use OTSU to threshold
rdna_thresh_array = multithresh(rdna_mat, 2); %cell and dark background removal
nuc_thresh_array = multithresh(nuc_mat, 2);

%% Apply threshold to image stacks
rdna_sig = rdna_mat;
rdna_sig(rdna_mat < rdna_thresh_array(end)) = nan;
nuc_sig = nuc_mat;
nuc_sig(nuc_mat < nuc_thresh_array(end)) = nan;

%% Statistical analysis
rho = corr(rdna_sig(:), nuc_sig(:), 'rows', 'complete');
rdna_vol = sum(~isnan(rdna_sig(:))) * pixel_size *pixel_size * z_step;
rdna_std = nanstd(rescale(rdna_sig(:)));
nuc_vol = sum(~isnan(nuc_sig(:))) * pixel_size * pixel_size * z_step;
nuc_std = nanstd(rescale(nuc_sig(:)));

%% Cdc14 Radial Analysis
%create binary images
rdna_bin = rdna_mat > rdna_thresh_array(end);
nuc_bin = nuc_mat > nuc_thresh_array(end);
rdna_T = regionprops3(rdna_bin, 'Volume', 'VoxelList');
nuc_T = regionprops3(nuc_bin, 'Volume', 'Centroid', 'SurfaceArea');
% filter tables to largest volume object only
[~, rdna_idx] = max(rdna_T.Volume);
[~, nuc_idx] = max(nuc_T.Volume);
rdna_coords = rdna_T.VoxelList{rdna_idx};
nuc_centroid = nuc_T.Centroid(nuc_idx, :);
subs = rdna_coords - repmat(nuc_centroid, [size(rdna_coords,1), 1]);
subs_nm = subs .* repmat([pixel_size, pixel_size, z_step], [size(subs, 1), 1]);
subs_nm_sq = subs_nm.^2;
displacement_nm = sqrt(sum(subs_nm_sq,2));
dispersion_ratio = var(displacement_nm)/mean(displacement_nm);
nuc_sph = pi()^(1/3)*(6*nuc_T.Volume(nuc_idx)).^(2/3)/nuc_T.SurfaceArea(nuc_idx);
