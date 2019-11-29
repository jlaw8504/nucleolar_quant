%%rnda_corr_summary Compares the Pearson correlation coefficients of the
%%lacO/LacI-GFP rDNA signal with CDC14-CFP and CBF5-RFP.
%
%NOTE : This summary script relies on relative locations of the directories
%so please do NOT move this script from current location

%% Define relative directory locations
s.cdc14.dir = fullfile(...
    'overlay', 'KBY6310p1_MH3342_SPC29RFP_CDC14CFP', 'all_max_proj');
s.cbf5.dir = fullfile(...
    'overlay', 'KBY6315p1_MH3342_SPC29RFP_CBF5RFP', 'all_max_proj');
s.smc4.dir = fullfile(...
    'overlay', 'KBY6392p1_CDC14GFP_SMC4mCherry', 'all_max_proj');
%% Define rdna_pattern and nuc_pattern variables
s.cdc14.rdna_pattern = '*GFP.tif';
s.cbf5.rdna_pattern = '*GFP.tif';
s.cdc14.nuc_pattern = '*CFP.tif';
s.cbf5.nuc_pattern = '*RFP.tif';
s.smc4.rdna_pattern = '*GFP.tif';
s.smc4.nuc_pattern = '*RFP.tif';
%% rdna_corr function calls
s.cdc14.corr_array = rdna_corr(s.cdc14.dir, s.cdc14.rdna_pattern, s.cdc14.nuc_pattern);
s.cbf5.corr_array = rdna_corr(s.cbf5.dir, s.cbf5.rdna_pattern, s.cbf5.nuc_pattern);
s.smc4.corr_array = rdna_corr(s.smc4.dir, s.smc4.rdna_pattern, s.smc4.nuc_pattern);
%% ttest2 comparison of correlation coefficients
[~, s.cdc14_vs_cbf5_pval] = ttest2(s.cdc14.corr_array, s.cbf5.corr_array);
[~, s.cdc14_vs_smc4_pval] = ttest2(s.cdc14.corr_array, s.smc4.corr_array);
%% Display two-color histograms
figure;
s.cdc14.hist = histogram(s.cdc14.corr_array, 'Normalization', 'probability');
hold on;
s.cbf5.hist = histogram(...
    s.cbf5.corr_array, 'Normalization', 'probability', 'BinWidth', ...
    s.cdc14.hist.BinWidth);
hold off;
legend('CDC14-CFP', 'CBF5-RFP');
xlabel('Pearson''s Correlation Coefficient');
ylabel('Frequency');