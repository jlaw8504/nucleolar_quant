load('nucleolus_frap_data.mat');
%only compare the exponents of Cbf5, Net1, and Rpa190
%MH3342 doesn't show recovery
%% filter out all signals with less than 20% recovery
S.cbf5.rec = S.cbf5.norm_mat(:,end) - S.cbf5.norm_mat(:,2);
S.net1.rec = S.net1.norm_mat(:,end) - S.net1.norm_mat(:,2);
S.rpa190.rec = S.rpa190.norm_mat(:,end) - S.rpa190.norm_mat(:,2);
% apply filter
S.cbf5.exponents(S.cbf5.rec < 0.2) = nan;
S.net1.exponents(S.net1.rec < 0.2) = nan;
S.rpa190.exponents(S.rpa190.rec < 0.2) = nan;

exps = [ ...
    S.cbf5.exponents; ...
    S.net1.exponents; ...
    S.rpa190.exponents ...
    ];
halflives = (log(2)./exps); %in seconds
group = [ ...
    repmat({'Cbf5-GFP'}, size(S.cbf5.exponents)); ...
    repmat({'Net1-GFP'}, size(S.net1.exponents)); ...
    repmat({'Rpa190-GFP'}, size(S.rpa190.exponents)) ...
    ];
[p, tbl, stats] = anova1(halflives, group);
%% Bar chart with error bars
means = [ ...
    nanmean(log(2)./S.cbf5.exponents), ...
    nanmean(log(2)./S.net1.exponents), ...
    nanmean(log(2)./S.rpa190.exponents) ...
    ];
sems = [ ...
    nanstd(log(2)./S.cbf5.exponents)/sqrt(numel(S.cbf5.exponents)), ...
    nanstd(log(2)./S.net1.exponents)/sqrt(numel(S.net1.exponents)), ...
    nanstd(log(2)./S.rpa190.exponents)/sqrt(numel(S.rpa190.exponents)) ...
    ];
figure;
bar(means);
hold on;
errorbar(means, sems, 'Marker', 'none', 'LineStyle', 'none', 'Color', 'k');
hold off;
xticklabels({'Cbf5-GFP', 'Net1-GFP', 'Rpa190-GFP'});