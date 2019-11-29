function s = tetra_dimer_summary
%%summary_quant_script Summarizes the cropped images of the Dimer and Tetra
%%lacO/LacI-GFP rDNA locus and CBF5.

%% Build summary structure array
root = pwd;
cd(root);

tags = {'LacO', 'CBF5'};
phase = {'G1', 'M'};
rdna_state = {'Dimer', 'Tetra'};
metrics = {'otsu_area_array', 'int_int_array', 'homog_array'};
for t_cnt = 1:numel(tags)
    s.(tags{t_cnt}) = struct;
    for p_cnt = 1:numel(phase)
        s.(tags{t_cnt}).(phase{p_cnt}) = struct;
        for r_cnt = 1:numel(rdna_state)
            s.(tags{t_cnt}).(phase{p_cnt}).(rdna_state{r_cnt}) = struct;
            dir_s = dir(sprintf('*_%s_%s-*FP*%s',...
                rdna_state{r_cnt}, tags{t_cnt}, phase{p_cnt}));
            if numel(dir_s) > 1 || numel(dir_s) == 0
                error('Found multiple or no directory matches');
            end
            [...
                s.(tags{t_cnt}).(phase{p_cnt}).(rdna_state{r_cnt}).(metrics{1}),...
                s.(tags{t_cnt}).(phase{p_cnt}).(rdna_state{r_cnt}).(metrics{2}),...
                s.(tags{t_cnt}).(phase{p_cnt}).(rdna_state{r_cnt}).(metrics{3}),...
                ] = quant_int(dir_s(1).name);
        end
    end
end

%% Statistical analysis portion
for t_cnt = 1:numel(tags)
    for p_cnt = 1:numel(phase)
        for m_cnt = 1:numel(metrics)
            [~, s.(tags{t_cnt}).(phase{p_cnt}).pvals.(metrics{m_cnt})] = ...
                ttest2(...
                s.(tags{t_cnt}).(phase{p_cnt}).Dimer.(metrics{m_cnt}),...
                s.(tags{t_cnt}).(phase{p_cnt}).Tetra.(metrics{m_cnt})...
                );
        end
    end
end

%% otsu_area_array parsing for JSON export
%Create individual LacO g and m structure arrays
g.dimer = s.LacO.G1.Dimer.otsu_area_array;
m.dimer = s.LacO.M.Dimer.otsu_area_array;
g.tetra = s.LacO.G1.Tetra.otsu_area_array;
m.tetra = s.LacO.M.Tetra.otsu_area_array;
fidG = fopen('quant_lacO_otsu2_data_G1.json', 'w');
fidM = fopen('quant_lacO_otsu2_data_M.json', 'w');
fprintf(fidG, '%s', jsonencode(g));
fprintf(fidM, '%s', jsonencode(m));
fclose(fidG);
fclose(fidM);

%Create individual CBF5 g and m structure arrays
g.dimer = s.CBF5.G1.Dimer.otsu_area_array;
m.dimer = s.CBF5.M.Dimer.otsu_area_array;
g.tetra = s.CBF5.G1.Tetra.otsu_area_array;
m.tetra = s.CBF5.M.Tetra.otsu_area_array;
fidG = fopen('quant_cbf5_otsu2_data_G1.json', 'w');
fidM = fopen('quant_cbf5_otsu2_data_M.json', 'w');
fprintf(fidG, '%s', jsonencode(g));
fprintf(fidM, '%s', jsonencode(m));
fclose(fidG);
fclose(fidM);