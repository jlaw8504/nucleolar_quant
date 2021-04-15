function [exponents, norm_mat, frap_mat] = frap_analysis_dir_mean_bleach(dir_name, num_timepoints)
%%frap_analysis_dir performs FRAP analysis on all ND2 files in a directory.
%
%   inputs :
%       dir : A string variable specifying the directory containing the ND2
%       sequence files to parse
%
%       num_timepoints : A scalar variable specifying the number of time
%       points in the timelapse following laser bleaching
%
%   outputs :
%       exponenets : A vector containing the exponential coefficients for
%       each normalized recovery curve, excluding the pre-bleach value.
%
%       norm_mat : A 2D matrix containing the normalized recovery curves.
%       Each row is a recovery curve and each column represents a time
%       point. The first value in the recovery curve, corresponding to the
%       pre-bleach intensity value, is set to 1.
%
%       frap_mat : A 2D matrix containing the background-subtracted,
%       photo-bleach corrected, mean intensity values for the bleached
%       region. Only the post-bleach values are included in this curve.
nd2_files =  dir(fullfile(dir_name, '*.nd2'));
% bleach_area_array = zeros([numel(nd2_files)/2, 1]);
frap_mat = zeros([numel(nd2_files)/2, num_timepoints]);
norm_mat = zeros([numel(nd2_files)/2, num_timepoints + 1]);
exponents = zeros([numel(nd2_files)/2, 1]);
idx = 1;
times = transpose([0:49]/30); %s
%% Calc mean bleach rate
bleach_slopes = zeros([numel(nd2_files)/2, 1]);
b_idx = 1;
for n = 1:2:numel(nd2_files)
    cell_seq = {nd2_files(n).name, nd2_files(n+1).name};
    [~, ~, post_stack] = parse_seq(cell_seq);
    bleach_slopes(b_idx) = bleach_corr_coeff(post_stack);
    b_idx = b_idx + 1;
end
corr_slope = mean(bleach_slopes(bleach_slopes < 0));
for n = 1:2:numel(nd2_files)
    cell_seq = {nd2_files(n).name, nd2_files(n+1).name};
    [pre_im, laser_im, post_stack] = parse_seq(cell_seq);
    [fg_bin, ~, ~, laser_bin, pre_val] = ...
        binary_process(pre_im, laser_im);
    if isnan(fg_bin)
        exponents(idx) = nan;
        norm_mat(idx,:) = nan;
        frap_mat(idx,:) = nan;
    else
%         area_bleach = 1 - sum(bg_bin(:))/sum(fg_bin(:));
        frap_means = simple_corr(post_stack, laser_bin, corr_slope);
        norm_means = frap_means/pre_val;
        exp_fit = fit(times, norm_means, 'exp1');
        coeffs = coeffvalues(exp_fit);
        exponents(idx) = coeffs(end);
        frap_mat(idx,:) = frap_means';
        norm_mat(idx,:) = [1, norm_means'];
%         bleach_area_array(idx) = area_bleach;
    end
    idx = idx + 1;
end