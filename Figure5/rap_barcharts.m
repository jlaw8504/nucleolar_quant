%% Calculate compression ratio (mean int/ vol)
load('rap_vol_ints.mat');
fnames = fieldnames(S);
treatments = {'dmso', 'rap'};
phases = {'g1', 'm'};
for f = 1:numel(fnames)
    for t = 1:numel(treatments)
        for p = 1:numel(phases)
            S.(fnames{f}).(treatments{t}).(phases{p}).comps = ...
                S.(fnames{f}).(treatments{t}).(phases{p}).mean_ints./ ...
                (S.(fnames{f}).(treatments{t}).(phases{p}).vols*(0.0645*0.0645* 0.2));
        end
    end
end
for f = 1:numel(fnames)
    S.(fnames{f}).g1_comp_p = ranksum( ...
        S.(fnames{f}).dmso.g1.comps, ...
        S.(fnames{f}).rap.g1.comps);
    S.(fnames{f}).m_comp_p = ranksum( ...
        S.(fnames{f}).dmso.m.comps, ...
        S.(fnames{f}).rap.m.comps);
end

[g1, m, labels] = parse_barchart_data(S);
%% G1
figure
bar(g1.mean_comps);
hold on;
errorbar( ...
    g1.mean_comps, ...
    g1.sem_comps, ...
    'Color', 'black', ...
    'LineStyle', 'none', ...
    'Marker', 'none');
hold off;
set(gca, 'XTickLabels', labels);
set(gca, 'XTickLabelRotation', 30);
title('G1');
ylabel('Mean Intensity/\mum^3');
%% Metaphase
figure
bar(m.mean_comps);
hold on;
errorbar( ...
    m.mean_comps, ...
    m.sem_comps, ...
    'Color', 'black', ...
    'LineStyle', 'none', ...
    'Marker', 'none');
hold off;
set(gca, 'XTickLabels', labels);
set(gca, 'XTickLabelRotation', 30);
title('Metaphase');
ylabel('Mean Intensity/\mum^3');