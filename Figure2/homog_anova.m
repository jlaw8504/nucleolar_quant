load homog_mip_sip.mat;
fnames = fieldnames(S);
mip_data_g1 = [];
mip_data_m = [];
group_g1 = {};
group_m = {};
sip_data_g1 = [];
sip_data_m = [];
for f = 1:numel(fnames)
    mip_data_g1 = [mip_data_g1; S.(fnames{f}).g1.mip_homog_array'];
    mip_data_m = [mip_data_m; S.(fnames{f}).m.mip_homog_array'];
    sip_data_g1 = [sip_data_g1; S.(fnames{f}).g1.sip_homog_array'];
    sip_data_m = [sip_data_m; S.(fnames{f}).m.sip_homog_array'];
    group_g1 = [group_g1; ...
        repmat(fnames(f), [numel(S.(fnames{f}).g1.sip_homog_array), 1])];
    group_m = [group_m; ...
        repmat(fnames(f), [numel(S.(fnames{f}).m.sip_homog_array), 1])];
end
sip_g1_c = my_chart(sip_data_g1, group_g1, fnames);
title('Sip G1');
sip_m_c = my_chart(sip_data_m, group_m, fnames);
title('Sip M');
mip_g1_c = my_chart(mip_data_g1, group_g1, fnames);
title('Mip G1');
mip_m_c = my_chart(mip_data_m, group_m, fnames);
title('Mip M');

function c= my_chart(data, group, fnames)
[~,~,stats] = anova1(data, group, 'off');
[c, m] = multcompare(stats, 'Display', 'off');
[means, sort_idx] = sort(m(:,1));
se = m(sort_idx, 2);
figure;
bar(means);
hold on;
errorbar(means, se, 'Color', 'black', 'LineStyle', 'none', 'Marker', 'none');
hold off;
set(gca, 'xticklabels', fnames(sort_idx));
set(gca, 'ylim', [0.5 0.75]);
ylabel('Mean Homogeneity Index');
end
