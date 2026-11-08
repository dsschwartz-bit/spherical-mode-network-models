%FILTER_EXCITATION_GAUSSIAN  Causality analysis of an S-parameter model and an excitation, with Gaussian filter.
%
%   1. S-parameter file: causality of every entry S_ij on its own
%      (check_snp_causality.m in matlab/snp): onset-lead grid, impulse
%      responses, worst entries, impulse vs filtered responses.
%   2. Excitation: the S-parameter block driven by the CST excitation,
%      y = IFFT(S X); non-causal energy of every significant response,
%      with all figures shown (excitationCausalityWorkflow.m).
%   3. Filter: if the responses miss the causality target (or with
%      force_filter = true), a delayed Gaussian is built, applied to the
%      excitation and re-checked, and the filtered signal is written for CST
%      (ASCII) and ADS (MDF) in OutDir.
%
% cst_to_ads_mdf.m only runs steps 2-3 when an excitation fails the target;
% use this script to look at the analysis or build a filter on demand.
% Inputs : CST time-signal export with the incident wave i<port> (in_file),
%          spherical-mode Touchstone model (snp_file, output/SnP).
% Outputs: figures; filtered CST/MDF files in output/cst_ads/filtered_excitations.

clear; clc;
here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', '..', '..', 'setup.m'));   % adds matlab/ (incl. matlab/snp) and third_party/
addpath(fullfile(here, 'utils')); % excitationCausalityWorkflow

%% ===== USER SETTINGS =====
cst_root = "C:\path\to\CST\";   % folder containing the CST project folders
% =========================

%% ---- settings -------------------------------------------------------------
snp_file = fullfile(repo_root, "output", "SnP", "S_dom_psi_rectangularpatch_symmetry_FDTD_CST.s14p");
in_file  = string(fullfile(cst_root, "rectangularpatch_symmetry_FDTD_2.5GHz", "Export", "time_signals", "port2_OC", "iowaves_gaussian_2GHz_3GHz.txt"));
% in_file = cst_root + "rectangularpatch_symmetry_FDTD\Export\time_signals\iowaves_gaussian_0.1GHz_6GHz.txt";   % alternative: broadband excitation
out_dir  = fullfile(repo_root, "output", "cst_ads", "filtered_excitations");
if ~isfolder(out_dir), mkdir(out_dir); end
port     = 1;                 % CST signal i<port> is the excitation

show_snp_analysis = true;     % step 1: causality figures of the S-parameter file itself
force_filter      = false;    % true: build and apply a filter even if the excitation meets the target
force_type        = "bandpass";  % forced filter: "bandpass" (Gaussian, -3 dB at force_band) or "lowpass"
force_band        = [];       % Hz [flo fhi] for the forced band-pass; [] = [lowest S frequency, fmax]

opts = { ...
    'EncTolPct',     0.01, ...   % % of a response's own energy allowed before t = 0
    'SigdB',         20, ...     % responses within this of the strongest are checked
    'FirstTargetdB', -40, ...    % first filter: excitation spectrum at fmax re its peak
    'StepdB',        10, ...     % if already below: this much lower per filter
    'ShiftSigmas',   5, ...      % Gaussian delay, in std. devs. of its impulse response
    'MaxIter',       6, ...
    'WriteIfCausal', false, ...  % also write (unfiltered) files when no filter is needed
    'PlotAlways',    true, ...   % show the causality figures even when the check passes
    'OutDir',        out_dir};

%% ---- 1. causality of the S-parameter file ------------------------------------
if show_snp_analysis
    runSnpCausality(snp_file, here);
end

%% ---- 2-3. excitation check, filter, CST / ADS files ----------------------------
out = excitationCausalityWorkflow(in_file, port, snp_file, opts{:}, ...
    'ForceFilter', force_filter, 'ForceType', force_type, 'ForceBand', force_band);

fprintf('\n');
if out.passed
    fprintf('Excitation meets the %g %% target (worst non-causal energy %.2g %%): no filter applied.\n', ...
        opts{find(strcmp(opts(1:2:end), 'EncTolPct'))*2}, out.encIn);
else
    if isempty(out.band), fd = sprintf('Gaussian low-pass, -3 dB at %.4g GHz', out.f3/1e9);
    else,                 fd = sprintf('Gaussian band-pass, -3 dB at %.4g and %.4g GHz', out.band/1e9); end
    fprintf('Filter: %s, delay %.1f ps (sigma_t %.1f ps)\n', fd, out.delay*1e12, out.sigma_t*1e12);
    fprintf('Worst non-causal energy: %.4g %% (input) -> %.4g %% (filtered)\n', out.encIn, out.encOut);
    if strlength(out.cstFile) > 0, fprintf('  CST: %s\n', out.cstFile); end
    if strlength(out.mdfFile) > 0, fprintf('  ADS: %s\n', out.mdfFile); end
end

%% ========================================================================
function runSnpCausality(snp_check_external, here) %#ok<INUSD> (snp_check_external is read by check_snp_causality)
% Run the script check_snp_causality.m (matlab/snp, on the path via setup.m)
% on the given Touchstone file in this function's workspace, so it does not
% clear or overwrite the variables of the calling script.
fprintf('S-parameter causality (check_snp_causality.m)\n');
check_snp_causality;
end
