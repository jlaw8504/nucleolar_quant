function [c, m] = jgl_anova1(data_cell, label_cell)

data = [];
group = {};
for n = 1:numel(data_cell)
    data = [data, data_cell{n}];
    labels = repmat(label_cell(n), size(data_cell{n}));
    group = [group, labels];
end
[~,~,stats] = anova1(data, group, 'off');
[c, m] = multcompare(stats, 'Display', 'off');