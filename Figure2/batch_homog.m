load('homog_dirs.mat', 'S');
fnames = fieldnames(S);
phases = {'g1', 'm'};
for f = 1:numel(fnames)
    for p = 1:numel(phases)
        [ ...
            S.(fnames{f}).(phases{p}).mip_homog_array, ...
            S.(fnames{f}).(phases{p}).sip_homog_array ...
            ] = quant_homog( ...
            S.(fnames{f}).(phases{p}).dir);
    end
end
clearvars -except S;
save('homog_mip_sip.mat');
clear S;