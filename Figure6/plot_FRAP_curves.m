%% Load in the data
S = load('nucleolus_frap_data.mat');
S = S.S;
%% Iterate over each fieldname
fnames = fieldnames(S);
for f = 1:numel(fnames)
    plot((30/50)*(0:50), mean(S.(fnames{f}).norm_mat)');
    if f == 1
        hold on;
    elseif f == numel(fnames)
        hold off;
    end
end
xlabel('Time (s)');
ylabel('Normalized Intensity');
legend({'rDNA-GFP', 'Cdc14-GFP', 'Net1-GFP', 'Rpa190-GFP', 'Cbf5-GFP'});