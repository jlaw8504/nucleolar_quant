load('rap_dirs.mat');
fnames = fieldnames(S);
treatments = {'dmso', 'rap'};
phases = {'g1', 'm'};
labels = cell([1,16]);
cnt = 1;
return_dir = pwd;
for f = 1:numel(fnames)
    for t = 1:numel(treatments)
        for p = 1:numel(phases)
            cd(return_dir);
            cd(S.(fnames{f}).(treatments{t}).(phases{p}).dir);
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
            S.(fnames{f}).(treatments{t}).(phases{p}).infocus_areas = outlier_free(:,1);
            S.(fnames{f}).(treatments{t}).(phases{p}).mip_areas = outlier_free(:,2);
            S.(fnames{f}).(treatments{t}).(phases{p}).sip_areas = outlier_free(:,3);
            S.(fnames{f}).(treatments{t}).(phases{p}).vols = outlier_free(:,4);
            S.(fnames{f}).(treatments{t}).(phases{p}).int_ints = outlier_free(:,5);
            S.(fnames{f}).(treatments{t}).(phases{p}).mean_ints = outlier_free(:,6);
            [ ...
                S.(fnames{f}).(treatments{t}).(phases{p}).vol_rho, ...
                S.(fnames{f}).(treatments{t}).(phases{p}).vol_p] = corr( ...
                S.(fnames{f}).(treatments{t}).(phases{p}).vols, ...
                S.(fnames{f}).(treatments{t}).(phases{p}).mean_ints, ...
                'type', 'Pearson', ...
                'rows', 'complete' ...
                );
            labels{cnt} = sprintf('%s %s %s', treatments{t}, phases{p}, fnames{f});
            vol_rhos(cnt) = S.(fnames{f}).(treatments{t}).(phases{p}).vol_rho;
            vol_ps(cnt) = S.(fnames{f}).(treatments{t}).(phases{p}).vol_p;
            clear outlier_mat;
            cnt = cnt + 1
        end
    end
end
cd(return_dir);