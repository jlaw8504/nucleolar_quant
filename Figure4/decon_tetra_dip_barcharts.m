%% Calculate compression ratio (mean int/ vol)
load('decon_dip_tet_data.mat');
fnames = fieldnames(S);
treatments = {'dimer', 'tetra'};
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
        S.(fnames{f}).dimer.g1.comps, ...
        S.(fnames{f}).tetra.g1.comps);
    S.(fnames{f}).m_comp_p = ranksum( ...
        S.(fnames{f}).dimer.m.comps, ...
        S.(fnames{f}).tetra.m.comps);
end

[g1, m, labels] = parse_dip_tetra_barchart_data(S);
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
title('Deconvolved, G1');
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
title('Deconvolved, Metaphase');
ylabel('Mean Intensity/\mum^3');