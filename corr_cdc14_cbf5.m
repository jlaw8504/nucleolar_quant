function [rho_array, cbf5_std_array, cbf5_area_array] = corr_cdc14_cbf5(cdc14_mip, cbf5_mip)
%% corr_cdc14_cbf5 calculates the correlation between two signals.
%
%
% inputs:
%  cdc14_mip : maximum intensity projection of the cdc14 signal.
%  
%  cbf5_mip : maximum intensity projection of the cbf5 signal.
%
% output:
%   rho : the correlation coefficient of both signals calculated between
%   each pair of columns.


%creates binary mask for CBF5 mip
thresh = multithresh(cbf5_mip,2);
binary = cbf5_mip>thresh(end);

%use imopen to eliminate smaller errorneous signals
binary2 = imopen(binary, strel('disk', 5));
bw_struct = bwconncomp(binary2);

%pre-allocate the rho_array and cbf5_std_array for this image
rho_array = zeros([bw_struct.NumObjects, 1]);
cbf5_std_array = rho_array;
cbf5_area_array = rho_array;

%Iterate over the bw_struct to generat multiple rhos
for n = 1:bw_struct.NumObjects
    cdc14_array = cdc14_mip(bw_struct.PixelIdxList{n});
    cbf5_array = cbf5_mip(bw_struct.PixelIdxList{n});
    cbf5_area_array(n,1) = numel(bw_struct.PixelIdxList{n});
    cbf5_std_array(n,1) = std(cbf5_array);
    rho_array(n,1) = corr(cdc14_array, cbf5_array);
end


