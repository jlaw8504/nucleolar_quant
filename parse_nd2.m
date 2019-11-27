function [cdc14_mip, cbf5_mip] = parse_nd2(nd2_filename, cdc14_ch, cbf5_ch)
%%parse_nd2 Parses an ND2 file and creates maximum intensity projections
%%for the cdc14 and cbf5 image channels.
%
%   inputs :
%       nd2_filename : A string variable pointing to the ND2 file you wish
%       to parse.
%
%       cdc14_ch : A scalar variable indicating the CDC14 channel number,
%       i.e 2.
%
%       cbf5_ch : A scalar variable indicating the CBF5 channel number,
%       i.e. 3.
%
%   outputs :
%       cdc14_mip : A 2D matrix containing the maximum intensity projection
%       of the CDC14 image channel.
%
%       cbf5_mip : A 2D matrix containing the maximum intensity projeciton
%       of the CBF5 image channel.

%% Sse bfopen to parse the nd2 file
im_cell = bfopen(nd2_filename);
%% I already know the dimensions of the image are
z_planes = 7;
%% pre-allocate the channel matrices
cdc14_im = zeros([size(im_cell{1,1}{1,1}), z_planes]);
cbf5_im = cdc14_im;
% instantiate channel counters
cdc14_cnt = 1;
cbf5_cnt = 1;
%iterate over im_cell{1,1} to parse images by channel
for n = 1:size(im_cell{1,1})
    if ~isempty(regexp(im_cell{1,1}{n,2}, ...
            sprintf('C=%d/3', cdc14_ch), 'once'))
        cdc14_im(:,:,cdc14_cnt) = im_cell{1,1}{n,1};
        cdc14_cnt = cdc14_cnt + 1;
    elseif ~isempty(regexp(im_cell{1,1}{n,2},...
            sprintf('C=%d/3', cbf5_ch), 'once'))
        cbf5_im(:,:,cbf5_cnt) = im_cell{1,1}{n,1};
        cbf5_cnt = cbf5_cnt + 1;
    end
end

% create max intensity projections for RFP and GFP channels
cdc14_mip = max(cdc14_im, [], 3);
cbf5_mip = max(cbf5_im, [], 3);


