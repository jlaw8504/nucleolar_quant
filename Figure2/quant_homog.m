function [mip_array, sip_array] = quant_homog(directory)
%%quant_int Use multithresh to separate signal from noise for CDC14 and
%%CBF5 fluorescent signals and perform quantitative intensity analysis on
%%the signal-labeled pixels.
%
%   input :
%       directory : A string variable pointing to the directory containing
%       the cropped, fluorescent images
%
%   outputs :
%       mip_array = An array variable containing the homogeniety values
%       generated from maximum intensity projections of the 3D image stacks
%
%       sip_array = An array variable containing the homogeniety values
%       generated from sum intensity projections of the 3D image stacks
%% Read in image
gfp_files = dir(fullfile(directory, '*GFP.tif'));
%% Pre-allocate array variables
mip_array = zeros([1,size(gfp_files,1)]);
sip_array = zeros([1,size(gfp_files,1)]);
for n = 1:size(gfp_files,1)
    filename = fullfile(gfp_files(n).folder,gfp_files(n).name);
    %open image with readTiffStack
    im = readTiffStack(filename);
    %convert image to double
    imdbl = double(im);
    %% Generate MIP and SIP
    %create a max-projection of the image (mip)
    mip = max(imdbl,[], 3);
    sip = sum(imdbl, 3);
    %% Find and NaN out zero values from MetaMorph Cropping
    if sum(mip(:) == 0) > 0 || sum(sip(:) == 0) > 0
        mip(mip == 0) = nan;
        sip(sip == 0) = nan;
        warning('%s contains values of zeros', filename);
    end
    %% Fix border problems
    %some borders were accidentally drawn on the margins of images
    %need to filter these out
    %get the locations of brightest pixels
    [~, idx] = max(mip(:));
    [row, col] = ind2sub(size(mip), idx);
    if row == 1 || col == 1 || row == size(mip,1) || col == size(mip, 2)
        warning('%s brightest pixel on border, excluding from analysis',...
            filename);
        mip_array(n) = nan;
        continue;
    end
    
    [~, idx] = max(sip(:));
    [row, col] = ind2sub(size(sip), idx);
    if row == 1 || col == 1 || row == size(sip,1) || col == size(sip, 2)
        warning('%s brightest pixel on border, excluding from analysis',...
            filename);
        sip_array(n) = nan;
        continue;
    end
    %% Background subtract and scale for mip
    thresh_mip = multithresh(mip,2);
    mip_rescale = mip;
    mip_otsu_bin = mip_rescale >= thresh_mip(end);
    mip_rescale(~mip_otsu_bin) = nan;
    mip_rescale = mip_rescale - min(mip_rescale(:));
    mip_rescale = mip_rescale./max(mip_rescale(:));
    %% Background subtract and scale for sip
    thresh_sip = multithresh(sip,2);
    sip_otsu_bin = sip >= thresh_sip(end);
    sip_rescale = sip;
    sip_rescale(~sip_otsu_bin) = nan;
    sip_rescale = sip_rescale - min(sip_rescale(:));
    sip_rescale = sip_rescale./max(sip_rescale(:));
    %% Calculate signal homogeneity using graycoprops for mip
    mip_gc_props = graycoprops(graycomatrix(mip_rescale), 'Homogeneity');
    sip_gc_props = graycoprops(graycomatrix(sip_rescale), 'Homogeneity');
    mip_array(n) = mip_gc_props.Homogeneity;
    sip_array(n) = sip_gc_props.Homogeneity;
end