%RUN_CAUSALITY_PORT2_OC_2P5GHZ  Causality check of one CST excitation against the S-parameter model.
% Runs excitationCausalityWorkflow for the 2.5 GHz patch with port 2 open
% (Gaussian 2-3 GHz excitation) against the spherical-mode Touchstone model,
% always showing the causality figures. No files are written unless a filter
% is needed (then into output/cst_ads/filtered_excitations).
% Inputs: CST incident-wave export (in_file) and the .s14p model in output/SnP.
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', '..', '..', 'setup.m'));   % adds matlab/ and third_party/
addpath(fullfile(here, 'utils')); % excitationCausalityWorkflow

%% ===== USER SETTINGS =====
in_file  = "C:\path\to\CST\rectangularpatch_symmetry_FDTD_2.5GHz\Export\time_signals\port2_OC\iowaves_gaussian_2GHz_3GHz.txt";
% =========================
snp_file = fullfile(repo_root,"output","SnP","S_dom_psi_rectangularpatch_symmetry_FDTD_CST.s14p");
out_dir  = fullfile(repo_root,"output","cst_ads","filtered_excitations");
if ~isfolder(out_dir), mkdir(out_dir); end
opts = {'EncTolPct',0.01,'SigdB',20,'FirstTargetdB',-40,'StepdB',10, ...
        'ShiftSigmas',5,'MaxIter',6,'WriteIfCausal',false,'PlotAlways',true, ...
        'OutDir',out_dir};
out = excitationCausalityWorkflow(in_file, 1, snp_file, opts{:});
disp(out)
