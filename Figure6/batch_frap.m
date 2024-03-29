%This script iteratively runs frap_analysis_dir function on hard-coded
%directories to generate a single structural array S
%% Hardcoded variables
image_file_root = 'C:\Users\lawrimor\University of North Carolina at Chapel Hill\Biology Bloom Lab - Documents\rDNA paper\FRAP';
folders = { ...
    'MH3342_FRAP_Good', ...
    'JLY1048_Cdc14GFP_good_frap', ...
    'net1_FRAP_good', ...
    'rpa190_FRAP_good', ...
    'cbf5_FRAP_good' ...
    };
dirs = cellfun(@(x) strcat(image_file_root, filesep, x),folders, 'UniformOutput', false); 
labels = { ...
    'rDNA-lacO/LacI-GFP', ...
    'Cdc14-GFP', ...
    'Net1-GFP', ...
    'Rpa190-GFP', ...
    'Cbf5-GFP' ...
    };
fnames = { ...
    'rdna', ...
    'cdc14', ...
    'net1', ...
    'rpa190', ...
    'cbf5' ...
    };
num_timepoints = 50;
for n = 1:numel(dirs)
    cd(dirs{n});
    [ ...
        S.(fnames{n}).exponents, ...
        S.(fnames{n}).norm_mat, ...
        S.(fnames{n}).frap_mat ...
        ] = frap_analysis_dir(pwd, num_timepoints);
    cd('..');
end
%% Anova1 analysis
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
[c, m] = multcompare(stats);
%% Bar chart with errorbars
figure;
bar(m(:,1));
hold on;
errorbar(m(:,1), m(:,2), 'LineStyle', 'none', 'Color', 'black');
hold off;

    