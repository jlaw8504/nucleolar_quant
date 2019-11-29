s = mutant_analysis(pwd);
%% otsu_area_array parsing for JSON export
%Create individual LacO g and m structure arrays
g.wt = s.CDC14.G1.wt.otsu_area_array;
m.wt = s.CDC14.M.wt.otsu_area_array;
g.ycg1d2 = s.CDC14.G1.ycg1d2.otsu_area_array;
m.ycg1d2 = s.CDC14.M.ycg1d2.otsu_area_array;
g.fob1d = s.CDC14.G1.fob1d.otsu_area_array;
m.fob1d = s.CDC14.M.fob1d.otsu_area_array;
fidG = fopen('json_CDC14_otsu2_data_G1.json', 'w');
fidM = fopen('json_CDC14_otsu2_data_M.json', 'w');
fprintf(fidG, '%s', jsonencode(g));
fprintf(fidM, '%s', jsonencode(m));
fclose(fidG);
fclose(fidM);

%Create individual CBF5 g and m structure arrays
g.wt = s.CBF5.G1.wt.otsu_area_array;
m.wt = s.CBF5.M.wt.otsu_area_array;
g.ycg1d2 = s.CBF5.G1.ycg1d2.otsu_area_array;
m.ycg1d2 = s.CBF5.M.ycg1d2.otsu_area_array;
g.fob1d = s.CBF5.G1.fob1d.otsu_area_array;
m.fob1d = s.CBF5.M.fob1d.otsu_area_array;
fidG = fopen('json_CBF5_otsu2_data_G1.json', 'w');
fidM = fopen('json_CBF5_otsu2_data_M.json', 'w');
fprintf(fidG, '%s', jsonencode(g));
fprintf(fidM, '%s', jsonencode(m));
fclose(fidG);
fclose(fidM);

%% ANOVA1 setup
tags = {'CDC14', 'CBF5'};
phase = {'G1', 'M'};
rdna_state = {'wt', 'ycg1d2','fob1d'};
metrics = {'otsu_area_array', 'int_int_array', 'log_array'};
for t_cnt = 1:numel(tags)
    s_anova.(tags{t_cnt}) = struct;
    for p_cnt = 1:numel(phase)
        s_anova.(tags{t_cnt}).(phase{p_cnt}) = struct;
        for m_cnt = 1:numel(metrics)
            s_anova.(tags{t_cnt}).(phase{p_cnt}).(metrics{m_cnt}) = struct;
            data = [];
            group = {};
            for r_cnt = 1:numel(rdna_state)
                data = [data;...
                    s.(tags{t_cnt}).(phase{p_cnt}).(rdna_state{r_cnt}).(metrics{m_cnt})'...
                    ];
                labels = repmat(rdna_state(r_cnt),...
                    [...
                    size(s.(tags{t_cnt}).(phase{p_cnt}).(rdna_state{r_cnt}).(metrics{m_cnt}), 2),1]);
                group = [group;labels];
            end
            s_anova.(tags{t_cnt}).(phase{p_cnt}).(metrics{m_cnt}).data = data;
            s_anova.(tags{t_cnt}).(phase{p_cnt}).(metrics{m_cnt}).group = group;
            [p, tbl, stats] = anova1(data, group, 'off');
            [c, m] = multcompare(stats, 'Display', 'off');
            s_anova.(tags{t_cnt}).(phase{p_cnt}).(metrics{m_cnt}).p = p;
            s_anova.(tags{t_cnt}).(phase{p_cnt}).(metrics{m_cnt}).tbl = tbl;
            s_anova.(tags{t_cnt}).(phase{p_cnt}).(metrics{m_cnt}).stats = stats;
            s_anova.(tags{t_cnt}).(phase{p_cnt}).(metrics{m_cnt}).c = c;
            s_anova.(tags{t_cnt}).(phase{p_cnt}).(metrics{m_cnt}).m = m;
            
        end
    end
end