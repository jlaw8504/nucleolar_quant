function s = quant_summary
%%quant_summary Build a structure array to contain all the intensity
%%measurements and plot the scatter and histograms comparing G1 and M phase

%% Create structure array
phases = {'G1' 'M'};
labels = {'CDC14', 'CBF5'};
for i = 1:numel(labels)
    for p = 1:numel(phases)
        s.(labels{i}).(phases{p}) = struct;
        %Find matching directories
        dir_string = sprintf('*_%s*_%s_55x55', labels{i}, phases{p});
        dir_struct = dir(dir_string);
        dir_name = dir_struct.name;
        %Gather data using quant_int function
        s.(labels{i}).(phases{p}).dir = dir_name;
        [s.(labels{i}).(phases{p}).otsu_area_array,...
            s.(labels{i}).(phases{p}).int_int_array,...
            s.(labels{i}).(phases{p}).homog_array] = ...
            quant_int(dir_name);
        %calculate mean intensity for each image
        s.(labels{i}).(phases{p}).mean_int_array = ...
            s.(labels{i}).(phases{p}).int_int_array...
            ./s.(labels{i}).(phases{p}).otsu_area_array;
        %% Scatter plots of mean int vs otsu area in nm
        figure;
        scatter(s.(labels{i}).(phases{p}).otsu_area_array * (0.0645^2),...
            s.(labels{i}).(phases{p}).mean_int_array);
        lsline;
        xlabel('Area (µm)');
        ylabel('Mean Intensity');
        title(sprintf('%s %s', labels{i}, phases{p}));
        %calculate rho and pval for corr between otsu area and mean
        %intensity values
        [s.(labels{i}).(phases{p}).area_mean_int_rho,...
            s.(labels{i}).(phases{p}).area_mean_int_pval] =...
            corr(...
            s.(labels{i}).(phases{p}).mean_int_array',...
            s.(labels{i}).(phases{p}).otsu_area_array');
        %calculate pval for ttest2 for homog_array between CDC14 and CBF5 for
        %each phase
        if strcmp('CBF5', labels{i})
            [~, s.compare.(phases{p}).pval] =...
                ttest2(...
                s.CDC14.(phases{p}).homog_array,...
                s.CBF5.(phases{p}).homog_array);
            %% Plot Histograms comparing sum of log correlation
            figure;
            cdc14_homog_hist =histogram(...
                s.CDC14.(phases{p}).homog_array, 'Normalization', 'probability');
            hold on;
            histogram(...
                s.CBF5.(phases{p}).homog_array, 'Normalization', 'probability',...
                'BinWidth', cdc14_homog_hist.BinWidth);
            hold off;
            legend({'CDC14-GFP', 'CBF5-GFP'}, 'Location', 'northwest');
            title(phases{p});
            xlabel('Signal Homogeneity');
            ylabel('Frequency');
        end
    end
end