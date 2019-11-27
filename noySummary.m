%% Build the summary structure
[s.wt.all_rhos, s.wt.all_stds, s.wt.all_cbf5_areas] = dir_rhos('KBY6317 ND', '*.nd2');
[s.gal.all_rhos, s.gal.all_stds, s.gal.all_cbf5_areas] = dir_rhos('KBY6413_NOY988 ND', 'gal*.nd2');
[s.glu.all_rhos, s.glu.all_stds, s.glu.all_cbf5_areas] = dir_rhos('KBY6413_NOY988 ND', '5hglu*.nd2');

%% Output the summary structure s to json file format for Python parsing
fid = fopen('cbf5.json', 'w');
fprintf(fid, '%s', jsonencode(s));
fclose(fid);

%% Generate ANOVA1 and multicomparison test
%Parse conditions and metrics
conditions = fieldnames(s);
metrics = fieldnames(s.(conditions{1}));
cnt_array = zeros([numel(conditions), 1]);
anova_p_array = zeros([numel(metrics), 1]);
multcompare_array = cell(size(anova_p_array));
for c = 1:numel(conditions)
    cnt_array(c,1) = numel(s.(conditions{c}).(metrics{1}));
end
for m = 1:numel(metrics)
    data = zeros([sum(cnt_array), 1]);
    group = cell([sum(cnt_array), 1]);
    start_idx = 1;
    for c = 1:numel(conditions)
        end_idx = sum(cnt_array(1:c));
        data(start_idx:end_idx) = s.(conditions{c}).(metrics{m});
        group(start_idx:end_idx) = repmat(conditions(c), [cnt_array(c),1]);
        start_idx = end_idx +1;
    end
    [anova_p_array(m), ~, stats] = anova1(data, group, 'off');
    multcompare_array{m} = multcompare(stats, 'Display', 'off');
end
        
            
        
%% Plot the comparison histograms
figure(1);
hist= histogram(s.glu.all_rhos, 'Normalization','probability');hold on;
histogram(s.gal.all_rhos, 'Normalization', 'probability','BinWidth', hist.BinWidth);
histogram(s.wt.all_rhos, 'Normalization', 'probability','BinWidth', hist.BinWidth);
hold off;
legend('glu','gal', 'wt');
xlabel('Pearson''s Correlation Coefficient CDC14 vs CBF5 Signal');
ylabel('Frequency');

figure(2)
hist= histogram(s.glu.all_cbf5_areas, 'Normalization','probability'); hold on
histogram(s.gal.all_cbf5_areas, 'Normalization', 'probability','BinWidth', hist.BinWidth);
histogram(s.wt.all_cbf5_areas, 'Normalization', 'probability','BinWidth', hist.BinWidth);
hold off;
legend('glu','gal','wt');
xlabel('Area of CBF5 Signal in Pixels');
ylabel('Frequency');

figure(3);
hist= histogram(s.glu.all_stds, 'Normalization','probability'); hold on
histogram(s.gal.all_stds, 'Normalization', 'probability','BinWidth', hist.BinWidth);
histogram(s.wt.all_stds, 'Normalization', 'probability','BinWidth', hist.BinWidth);
hold off;
legend('glu','gal', 'wt');
xlabel('Standard Deviation of CBF5 Signal');
ylabel('Frequency');