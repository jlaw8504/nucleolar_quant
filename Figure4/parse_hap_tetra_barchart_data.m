function [g1, m, labels] = parse_hap_tetra_barchart_data(S)
%%parse_barchart_data parses structured array to gather mean ans SEM
%%statistics on DMSO and Rapamycin treated cells

labels = { ...
    'Nop56-RFP Dimer', ...
    'Nop56-RFP Tetramer', ...
    'Net1-RFP Dimer', ...
    'Net1-RFP Tetramer' ...
    };
%% Vols
g1.mean_vols = [ ...
    nanmean(S.nop56.dimer.g1.vols), ...
    nanmean(S.nop56.tetra.g1.vols), ...
    nanmean(S.net1.dimer.g1.vols), ...
    nanmean(S.net1.tetra.g1.vols), ...
    ];

g1.sem_vols = [ ...
    nanstd(S.nop56.dimer.g1.vols)/sqrt(numel(S.nop56.dimer.g1.vols)), ...
    nanstd(S.nop56.tetra.g1.vols)/sqrt(numel(S.nop56.tetra.g1.vols)), ...
    nanstd(S.net1.dimer.g1.vols)/sqrt(numel(S.net1.dimer.g1.vols)), ...
    nanstd(S.net1.tetra.g1.vols)/sqrt(numel(S.net1.tetra.g1.vols)) ...
    ];

m.mean_vols = [ ...
    nanmean(S.nop56.dimer.m.vols), ...
    nanmean(S.nop56.tetra.m.vols), ...
    nanmean(S.net1.dimer.m.vols), ...
    nanmean(S.net1.tetra.m.vols) ...
    ];

m.sem_vols = [ ...
    nanstd(S.nop56.dimer.m.vols)/sqrt(numel(S.nop56.dimer.m.vols)), ...
    nanstd(S.nop56.tetra.m.vols)/sqrt(numel(S.nop56.tetra.m.vols)), ...
    nanstd(S.net1.dimer.m.vols)/sqrt(numel(S.net1.dimer.m.vols)), ...
    nanstd(S.net1.tetra.m.vols)/sqrt(numel(S.net1.tetra.m.vols)) ...
    ];
%% Mean ints
g1.mean_ints = [ ...
    nanmean(S.nop56.dimer.g1.mean_ints), ...
    nanmean(S.nop56.tetra.g1.mean_ints), ...
    nanmean(S.net1.dimer.g1.mean_ints), ...
    nanmean(S.net1.tetra.g1.mean_ints), ...
    ];

g1.sem_mean_ints = [ ...
    nanstd(S.nop56.dimer.g1.mean_ints)/sqrt(numel(S.nop56.dimer.g1.mean_ints)), ...
    nanstd(S.nop56.tetra.g1.mean_ints)/sqrt(numel(S.nop56.tetra.g1.mean_ints)), ...
    nanstd(S.net1.dimer.g1.mean_ints)/sqrt(numel(S.net1.dimer.g1.mean_ints)), ...
    nanstd(S.net1.tetra.g1.mean_ints)/sqrt(numel(S.net1.tetra.g1.mean_ints)) ...
    ];

m.mean_ints = [ ...
    nanmean(S.nop56.dimer.m.mean_ints), ...
    nanmean(S.nop56.tetra.m.mean_ints), ...
    nanmean(S.net1.dimer.m.mean_ints), ...
    nanmean(S.net1.tetra.m.mean_ints) ...
    ];

m.sem_mean_ints = [ ...
    nanstd(S.nop56.dimer.m.mean_ints)/sqrt(numel(S.nop56.dimer.m.mean_ints)), ...
    nanstd(S.nop56.tetra.m.mean_ints)/sqrt(numel(S.nop56.tetra.m.mean_ints)), ...
    nanstd(S.net1.dimer.m.mean_ints)/sqrt(numel(S.net1.dimer.m.mean_ints)), ...
    nanstd(S.net1.tetra.m.mean_ints)/sqrt(numel(S.net1.tetra.m.mean_ints)) ...
    ];
%% Comps
g1.mean_comps =  [ ...
    nanmean(S.nop56.dimer.g1.comps), ...
    nanmean(S.nop56.tetra.g1.comps), ...
    nanmean(S.net1.dimer.g1.comps), ...
    nanmean(S.net1.tetra.g1.comps), ...
    ];

g1.sem_comps =  [ ...
    nanstd(S.nop56.dimer.g1.comps)/sqrt(numel(S.nop56.dimer.g1.comps)), ...
    nanstd(S.nop56.tetra.g1.comps)/sqrt(numel(S.nop56.tetra.g1.comps)), ...
    nanstd(S.net1.dimer.g1.comps)/sqrt(numel(S.net1.dimer.g1.comps)), ...
    nanstd(S.net1.tetra.g1.comps)/sqrt(numel(S.net1.tetra.g1.comps)) ...
    ];

m.mean_comps =  [ ...
    nanmean(S.nop56.dimer.m.comps), ...
    nanmean(S.nop56.tetra.m.comps), ...
    nanmean(S.net1.dimer.m.comps), ...
    nanmean(S.net1.tetra.m.comps) ...
    ];

m.sem_comps =  [ ...
    nanstd(S.nop56.dimer.m.comps)/sqrt(numel(S.nop56.dimer.m.comps)), ...
    nanstd(S.nop56.tetra.m.comps)/sqrt(numel(S.nop56.tetra.m.comps)), ...
    nanstd(S.net1.dimer.m.comps)/sqrt(numel(S.net1.dimer.m.comps)), ...
    nanstd(S.net1.tetra.m.comps)/sqrt(numel(S.net1.tetra.m.comps)) ...
    ];