function s = pNOY130_area_summary
%% Pattern parts
condition = {'GAL', 'GLU', 'WT'};
marker = {'Cbf5mCh', 'Cdc14GFP'};
phase = {'G1', 'M'};
%% Loop through parameters to parse images and validate sizes
for c = 1:numel(condition)
    for p = 1:numel(phase)
        for m = 1:numel(marker)
            s.(condition{c}).(phase{p}).(marker{m}).folder = ...
                sprintf('%s 55x55 %s %s', condition{c}, marker{m}, phase{p});
            s.(condition{c}).(phase{p}).(marker{m}).files = ...
                dir(fullfile(...
                s.(condition{c}).(phase{p}).(marker{m}).folder, '*.tif'));
            for i = 1:numel(s.(condition{c}).(phase{p}).(marker{m}).files)
                info = imfinfo(fullfile(...
                    s.(condition{c}).(phase{p}).(marker{m}).folder,...
                    s.(condition{c}).(phase{p}).(marker{m}).files(i).name));
                if numel(info) ~= 7 || info(1).Height ~= 55 || info(1).Width ~= 55
                    error('%s is wrong!', fullfile(...
                        s.(condition{c}).(phase{p}).(marker{m}).folder,...
                        s.(condition{c}).(phase{p}).(marker{m}).files(i).name));
                else
                    s.(condition{c}).(phase{p}).(marker{m}).stacks(:,:,:,i) = ...
                        readTiffStack( fullfile(...
                        s.(condition{c}).(phase{p}).(marker{m}).folder,...
                        s.(condition{c}).(phase{p}).(marker{m}).files(i).name));
                end
                
            end
            %% Create Max Intensity Projections
            s.(condition{c}).(phase{p}).(marker{m}).mips = ...
                squeeze(max(s.(condition{c}).(phase{p}).(marker{m}).stacks, [], 3));
            if m == 2
                if size(s.(condition{c}).(phase{p}).(marker{m}).mips,3) ~= ...
                        size(s.(condition{c}).(phase{p}).(marker{m-1}).mips,3)
                    error('%s 55x55 %s marker images different', condition{c}, phase{p});
                end
                for n = 1:size(s.(condition{c}).(phase{p}).(marker{m}).mips, 3)
                    mip_cbf5 = s.(condition{c}).(phase{p}).(marker{m-1}).mips(:,:, n);
                    mip_cdc14 = s.(condition{c}).(phase{p}).(marker{m}).mips(:,:, n);
                    % calculate cbf5 binary from Otsu threshold
                    cbf5_thresh = multithresh(mip_cbf5, 2);
                    cbf5_bin = mip_cbf5 > cbf5_thresh(end);
                    % calculate cdc14 binary from Otsu threshold
                    cdc14_thresh = multithresh(mip_cdc14, 2);
                    cdc14_bin = mip_cdc14 > cdc14_thresh(end);
                    % cbf5 to cdc14 signal correlation
                    s.(condition{c}).(phase{p}).corr(n) = ...
                        corr(mip_cbf5(cbf5_bin), mip_cdc14(cbf5_bin));
                    % store the area and std information
                    s.(condition{c}).(phase{p}).(marker{1}).areas(n) = sum(cbf5_bin(:));
                    s.(condition{c}).(phase{p}).(marker{2}).areas(n) = sum(cdc14_bin(:));
                    s.(condition{c}).(phase{p}).(marker{1}).stds(n) = std(mip_cbf5(cbf5_bin));
                    s.(condition{c}).(phase{p}).(marker{2}).stds(n) = std(mip_cdc14(cdc14_bin));
%                     if c == 3
%                         h = figure('WindowState', 'maximized');
%                         subplot(1,2,1);
%                         imshow(mip_cdc14, []);
%                         subplot(1,2,2);
%                         imshow(cdc14_bin);
%                         waitforbuttonpress;
%                         close(h);
%                     end
                end
            end
        end
    end
end

%% Write out JSON files
cdc14.wt = [s.WT.G1.Cdc14GFP.areas, s.WT.M.Cdc14GFP.areas]*(0.0645^2);
cdc14.wt = cdc14.wt(~isoutlier(cdc14.wt));
cdc14.gal = [s.GAL.G1.Cdc14GFP.areas, s.GAL.M.Cdc14GFP.areas]*(0.0645^2);
cdc14.gal = cdc14.gal(~isoutlier(cdc14.gal));
cdc14.glu = [s.GLU.G1.Cdc14GFP.areas, s.GLU.M.Cdc14GFP.areas]*(0.0645^2);
cdc14.glu = cdc14.glu(~isoutlier(cdc14.glu));
fid = fopen('json_CDC14_otsu2_data_G1.json', 'w');
fprintf(fid, '%s', jsonencode(cdc14));
fclose(fid);

cbf5.wt = [s.WT.G1.Cbf5mCh.areas, s.WT.M.Cbf5mCh.areas]*(0.0645^2);
cbf5.wt = cbf5.wt(~isoutlier(cbf5.wt));
cbf5.gal = [s.GAL.G1.Cbf5mCh.areas, s.GAL.M.Cbf5mCh.areas]*(0.0645^2);
cbf5.gal = cbf5.gal(~isoutlier(cbf5.gal));
cbf5.glu = [s.GLU.G1.Cbf5mCh.areas, s.GLU.M.Cbf5mCh.areas]*(0.0645^2);
cbf5.glu = cbf5.glu(~isoutlier(cbf5.glu));
fid = fopen('json_CBF5_otsu2_data_G1.json', 'w');
fprintf(fid, '%s', jsonencode(cbf5));
fclose(fid);
