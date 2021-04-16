%% Load in the data
load('nucleolus_frap_data.mat');
S = orderfields(S,  {'rdna', 'cdc14', 'net1', 'rpa190', 'cbf5'});
%% Iterate over each fieldname
fnames = fieldnames(S);
data_cell = cell([numel(fnames), 1]);
label_cell = cell([numel(fnames), 1]);
for f = 1:numel(fnames)
    norm_mat = S.(fnames{f}).norm_mat;
    perc_rec_array = norm_mat(:,end) - norm_mat(:,2);
    data_cell{f} = perc_rec_array * 100;
    S.(fnames{f}).perc_rec = perc_rec_array;
    label_cell{f} = repmat(fnames(f), size(perc_rec_array));
end
[p, tbl, stats] = anova1(vertcat(data_cell{:}), vertcat(label_cell{:}));
figure;
multcompare(stats)