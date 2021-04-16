load('homog_dirs.mat');
fnames = fieldnames(S);
phases = {'g1', 'm'};
vol_ps = zeros([numel(fnames), 1]);
vol_rhos = zeros([numel(fnames), 1]);
labels = cell([1,6]);
cnt = 1;
return_dir = pwd;
for f = 1:numel(fnames)
    for p = 1:numel(phases)
        display(cnt);
        cd(return_dir)
        cd(S.(fnames{f}).(phases{p}).dir);
        [ ...
            outlier_mat(:,1), ...
            outlier_mat(:,2), ...
            outlier_mat(:,3), ...
            outlier_mat(:,4), ...
            outlier_mat(:,5), ...
            outlier_mat(:,6) ...
            ] = quant_areas_vols('*GFP.tif');
        out_idx = sum(isoutlier(outlier_mat), 2) > 0;
        outlier_free = outlier_mat(~out_idx, :);
        S.(fnames{f}).(phases{p}).infocus_areas = outlier_free(:,1);
        S.(fnames{f}).(phases{p}).mip_areas = outlier_free(:,2);
        S.(fnames{f}).(phases{p}).sip_areas = outlier_free(:,3);
        S.(fnames{f}).(phases{p}).vols = outlier_free(:,4);
        S.(fnames{f}).(phases{p}).int_ints = outlier_free(:,5);
        S.(fnames{f}).(phases{p}).mean_ints = outlier_free(:,6);
        [ ...
            S.(fnames{f}).(phases{p}).vol_rho, ...
            S.(fnames{f}).(phases{p}).vol_p] = corr( ...
            S.(fnames{f}).(phases{p}).vols, ...
            S.(fnames{f}).(phases{p}).mean_ints, ...
            'type', 'Pearson', ...
            'rows', 'complete' ...
            );
        labels{cnt} = sprintf('%s %s', phases{p}, fnames{f});
        vol_rhos(cnt) = S.(fnames{f}).(phases{p}).vol_rho;
        vol_ps(cnt) = S.(fnames{f}).(phases{p}).vol_p;
        clear outlier_mat;
        cnt = cnt + 1;
    end
end
cd(return_dir);
% create Cdc14-GFP G1 and Cbf5-GFP, G1 mean intensity vs volume scatter
% plots
voxel_vol_um3 = (0.0645 * 0.0645 * 0.2);
figure;
scatter(S.cdc14.g1.vols * voxel_vol_um3, S.cdc14.g1.mean_ints);
xlabel('Volume(\mu^3)');
ylabel('Mean Intensity Above Background (AU)');
title('Cdc14-GFP G1')
figure;
scatter(S.cbf5.g1.vols * voxel_vol_um3, S.cbf5.g1.mean_ints);
xlabel('Volume(\mu^3)');
ylabel('Mean Intensity Above Background (AU)');
title('Cbf5-GFP G1')
% create bar charts
labels = {...
    'Cdc14-GFP', ...
    'Net1-GFP', ...
    'Gar1-GFP', ...
    'Nop56-GFP', ...
    'Nop1-GFP', ...
    'Cbf5-GFP'...
    };
% G1 barchart
figure;
bar(vol_rhos(1:2:end));
title('G1');
xticklabels(labels);
xtickangle(45)
ylabel('Pearson Correlation');
% G2/M barchart
figure;
bar(vol_rhos(2:2:end));
title('G2/M');
xticklabels(labels);
xtickangle(45);
ylabel('Pearson Correlation');