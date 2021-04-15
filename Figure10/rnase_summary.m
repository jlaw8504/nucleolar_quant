%% Cdc14-GFP Cbf5-mCherry file parsing
load rnaseH_dirs.mat;
treatments = {'vector', 'rnaseH'};
phases = {'g1', 'm'};
for t = 1:numel(treatments)
    for p = 1:numel(phases)
        cd(S.(treatments{t}).(phases{p}).dir);
        cdc14_files = dir('*GFP.tif');
        for f = 1:numel(cdc14_files)
            rdna_filename = cdc14_files(f).name;
            nuc_filename = strrep(rdna_filename, '_GFP.tif', '_RFP.tif');
            [ ...
                S.(treatments{t}).(phases{p}).rhos(f), ...
                S.(treatments{t}).(phases{p}).rdna_vols(f), ...
                S.(treatments{t}).(phases{p}).rdna_stds(f), ...
                S.(treatments{t}).(phases{p}).nuc_vols(f), ...
                S.(treatments{t}).(phases{p}).nuc_stds(f), ...
                S.(treatments{t}).(phases{p}).dispersion_ratio(f), ...
                S.(treatments{t}).(phases{p}).nuc_sph(f) ...
                ] = rnase_quant(rdna_filename, nuc_filename, 0.0645, 0.2);
        end
        %summarize rhos
        S.(treatments{t}).(phases{p}).mean_rhos = mean( ...
            S.(treatments{t}).(phases{p}).rhos);
        S.(treatments{t}).(phases{p}).sem_rhos = std( ...
            S.(treatments{t}).(phases{p}).rhos)./ ...
            sqrt(numel(S.(treatments{t}).(phases{p}).rhos));
        %summarize rdna vols
        S.(treatments{t}).(phases{p}).mean_rdna_vols = mean( ...
            S.(treatments{t}).(phases{p}).rdna_vols);
        S.(treatments{t}).(phases{p}).sem_rdna_vols = std( ...
            S.(treatments{t}).(phases{p}).rdna_vols)./ ...
            sqrt(numel(S.(treatments{t}).(phases{p}).rdna_vols));
        %sumarize rdna stds
        S.(treatments{t}).(phases{p}).mean_rdna_stds = mean( ...
            S.(treatments{t}).(phases{p}).rdna_stds);
        S.(treatments{t}).(phases{p}).sem_rdna_stds = std( ...
            S.(treatments{t}).(phases{p}).rdna_stds)./ ...
            sqrt(numel(S.(treatments{t}).(phases{p}).rdna_stds));
        %summarize nuc vols
        S.(treatments{t}).(phases{p}).mean_nuc_vols = mean( ...
            S.(treatments{t}).(phases{p}).nuc_vols);
        S.(treatments{t}).(phases{p}).sem_nuc_vols = std( ...
            S.(treatments{t}).(phases{p}).nuc_vols)./ ...
            sqrt(numel(S.(treatments{t}).(phases{p}).nuc_vols));
        %summarize nuc stds
        S.(treatments{t}).(phases{p}).mean_nuc_stds = mean( ...
            S.(treatments{t}).(phases{p}).nuc_stds);
        S.(treatments{t}).(phases{p}).sem_nuc_stds = std( ...
            S.(treatments{t}).(phases{p}).nuc_stds)./ ...
            sqrt(numel(S.(treatments{t}).(phases{p}).nuc_stds));
        
    end
end
%% Cdc14-GFP Cbf5-mCherry bar charts
metrics = {'rhos', 'nuc_vols', 'nuc_stds', 'rdna_vols', 'rdna_stds', 'dispersion_ratio', 'nuc_sph'};
for m = 1:numel(metrics)
    figure;
    [S.(metrics{m}).c_g1, m_g1] = jgl_anova1({...
        S.vector.g1.(metrics{m}), ...
        S.rnaseH.g1.(metrics{m}) ...
        }, ...
        {'G1, Vector', ...
        'G1, RnaseH' ...
        });
    [S.(metrics{m}).c_m, m_m] = jgl_anova1({...
        S.vector.m.(metrics{m}), ...
        S.rnaseH.m.(metrics{m}) ...
        }, ...
        {'G2/M, Vector', ...
        'G2/M, RnaseH' ...
        });
    if contains(metrics{m}, 'nuc')
        b = bar([m_g1(:,1); m_m(:,1)], 'LineWidth', 1);
        b.FaceColor = 'flat';
        b.CData = [repmat([178/255, 34/255, 34/255], [2,1]); repmat([128/255, 0, 0], [2,1])];
    elseif contains(metrics{m}, 'rdna')
        b = bar([m_g1(:,1); m_m(:,1)]);
        b.FaceColor = 'flat';
        b.CData = [repmat([154/255, 205/255, 50/255], [2, 1]); repmat([0/255, 128/255, 0/255], [2, 1])];
    else
        b = bar([m_g1(:,1); m_m(:,1)]);
        b.FaceColor = 'flat';
        b.CData = [repmat([0/255, 178/255, 178/255], [2,1]); repmat([0/255, 128/255, 128/255], [2, 1])];
    end
    hold on;
    errorbar([m_g1(:,1); m_m(:,1)], [m_g1(:,2); m_m(:,2)], ...
        'LineStyle', 'none', 'Marker', 'none', 'Color', 'black', 'LineWidth', 1);
    hold off
    xlim([0.5, 4.5]);
    switch m
        case 1
            ylabel({'Cdc14-GFP:Cbf5-mCherry','Correlation'});
        case 2
            ylabel('Cbf5-mCherry Volume (\mum^3');
        case 3
            ylabel('Rescaled Cbf5-mCherry Signal STD');
        case 4
            ylabel('Cdc14-GFP Volume (\mum^3)');
        case 5
            ylabel('Rescaled Cdc14-GFP Signal STD');
    end
    xticklabels(...
        {'G1, Vector', 'G1,RnaseH', ...
        'G2/M, Vector', 'G2/M, RnaseH'});
    set(gcf, 'Units', 'inches');
    if strcmp(metrics{m}, 'rhos')
        set(gcf, 'Position', [5,5,6.5, 2.5]);
        set(gca, 'XTickLabelRotation', 25);
        set(gca, 'LineWidth', 1)
    else
        set(gcf, 'Position', [5,5,3.25, 2.5]);
        set(gca, 'XTickLabelRotation', 25);
        set(gca, 'LineWidth', 1);
    end
    InSet = get(gca, 'TightInset');
    set(gca, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
    
end

