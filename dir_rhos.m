function [all_rhos, all_stds, all_cbf5_areas] = dir_rhos(directory, pattern)
%%dir_rhos Parses a directory for all image files that match a provided
%%pattern, generates a maximum intensity projection of channels 2 and 3,
%%generates a binary mask of channel 2 and applies it to channels 2 and 3,
%%and then calculates Pearson's correlation coefficient between channels 2
%%and 3 and stores them in a cell array.


%% Parse files in directory with pattern
files = dir(fullfile(directory, pattern));
%% Pre-allocate rho_array
rho_cell = cell([numel(files),1]);
std_cell = rho_cell;
area_cell = rho_cell;
%% Iterate over all files
for n = 1:numel(files)
    fullfilename = fullfile(files(n).folder, files(n).name);
    %assuming cdc14 channel is 2 and cbf5 channel is 3 (out of 3 channels)
    [cdc14_mip, cbf5_mip] = parse_nd2(fullfilename, 2, 3);
    [rho_cell{n}, std_cell{n}, area_cell{n}] = corr_cdc14_cbf5(cdc14_mip, cbf5_mip);
end
all_rhos = cell2mat(rho_cell);
all_stds = cell2mat(std_cell);
all_cbf5_areas = cell2mat(area_cell);