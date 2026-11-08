%GENERATE_RATIONALFIT_SPICE  Rational fit of a Touchstone file and SPICE netlist export.
% Fits the S-parameters of an SnP file (e.g. S_dom_*.sNp written by the
% spherical-mode workflows) with rational() or rationalfit(), plots data vs
% fit for the worst-fitting couplings, and writes:
%   <rfit_dir>/<filename>[_<PoleSharing>].ckt       SPICE netlist (generateSPICE),
%       node names renamed to the SnP port names and gmin resistors added;
%   <rfit_dir>/<filename>[_<PoleSharing>]_Rfit.mat  the fit itself (for Simulink).
% Requires RF Toolbox (rational() needs R2020a or later).

here = fileparts(mfilename('fullpath'));
run(fullfile(here, '..', '..', 'setup.m'));   % adds matlab/ and third_party/

%% ===== USER SETTINGS =====
filename="S_dom_Q_bowtie_edgeport_45degrees_PW_excitation.s6p";
% Alternatives:
%filename="S_dom_psi_rectangularpatch_symmetry_FDTD_CST.s14p";
%filename="S_dom_Q_loop_noloads.s10p";
filepath=string(fullfile(repo_root,'output','SnP'));                 % folder containing the SnP file
rfit_dir=string(fullfile(repo_root,'output','SnP','rationalfit'));   % netlist / fit output folder
% =========================
if ~exist(rfit_dir,'dir'), mkdir(rfit_dir); end

%% Fit settings
% Spherical mode ports to plot, by port number in the SnP file. Each one is
% plotted against every antenna (non-mode) port, e.g. [4 5] with antenna
% ports 1 and 2 plots S_{1,4}, S_{1,5}, S_{2,4} and S_{2,5}. Leave empty to
% plot only the worst-fit couplings (plus the antenna ports if
% always_plot_antenna_ports is true).
mode_ports_to_plot=[];

S=sparameters(fullfile(filepath,filename))

tolerance_dB=-40;
% Pole sharing: "All" = one pole set for every entry (fewest poles),
% "Column" = one set per column (excitation port), "Individual" = one set
% per entry (most poles). If the circuit simulation with "All" is
% inaccurate, try "Column" or "Individual".
pole_sharing="All";
% Error metric used to rank the couplings (and shown on the plots), in dB per entry:
% "absolute" = 10*log10(mean|fit-data|^2);
% "relative" = 10*log10(sum|fit-data|^2/sum|data|^2).
error_metric="absolute";
% Ports of the n_worst worst-fitting couplings (any port pair, including
% S11, S22) are always plotted.
n_worst=3;
% Also plot the antenna ports (S11, S22, S21, S12) regardless of fit error.
always_plot_antenna_ports=true;
% Fitting function: "rational" = rational() (RF Toolbox R2020a+, usually
% fewer poles for the same error), "rationalfit" = rationalfit(). If the
% rational() fit has a relative error above fallback_dB in any entry, the
% fit is redone with rationalfit().
fit_method="rational";
max_poles=100;                % pole cap (rational: MaxPoles; rationalfit: NPoles=[1 max_poles])
% Worst relative entry error, 10*log10(sum|fit-data|^2/sum|data|^2) (all-zero
% entries ignored), above which the rational() fit is rejected.
fallback_dB=tolerance_dB+10;  % e.g. -30 dB for a -40 dB tolerance

%% Rational fit
fit_rel_err=@(R) local_worst_rel_err(R,S);
switch fit_method
    case "rational"
        [Rfit,Rfit_errordB]=rational(S,...
            "Tolerance",tolerance_dB,"TendsToZero",true,...
            "PoleSharing",pole_sharing,"MaxPoles",max_poles);
        [worst_dB,worst_ij]=fit_rel_err(Rfit);
        fprintf('rational(), pole sharing "%s": worst entry S_{%d,%d}, relative error %.1f dB (fallback above %.1f dB).\n',...
            pole_sharing,worst_ij(1),worst_ij(2),worst_dB,fallback_dB);
        if worst_dB>fallback_dB
            warning('rational() fit inaccurate (S_{%d,%d}: %.1f dB > %.1f dB): refitting with rationalfit().',worst_ij(1),worst_ij(2),worst_dB,fallback_dB);
            fit_method="rationalfit";
        end
    case "rationalfit"
    otherwise
        error('fit_method must be "rational" or "rationalfit".');
end
if fit_method=="rationalfit"
    [Rfit,Rfit_errordB]=rationalfit(S,...
        "Tolerance",tolerance_dB,...
        "TendsToZero",true,"IterationLimit",[4 100],...
        "PoleSharing",pole_sharing,"NPoles",[1 max_poles]);
    [worst_dB,worst_ij]=fit_rel_err(Rfit);
    fprintf('rationalfit(), pole sharing "%s": worst entry S_{%d,%d}, relative error %.1f dB.\n',pole_sharing,worst_ij(1),worst_ij(2),worst_dB);
end
Rfit_errordB=max(Rfit_errordB(:));   % one number (rational() returns one per entry / column without shared poles)

% Number of distinct poles (entries sharing a pole set each carry a copy, so do not sum over all entries):
if isa(Rfit,'rational')
    P=Rfit.Poles;
    switch pole_sharing
        case "All",    n_poles=numel(P);
        case "Column", n_poles=sum(cellfun(@numel,P(1,:)));
        otherwise,     n_poles=sum(cellfun(@numel,P(:)));
    end
else
    switch pole_sharing
        case "All",    n_poles=numel(Rfit(1,1).A);
        case "Column", n_poles=sum(arrayfun(@(j) numel(Rfit(1,j).A),1:S.NumPorts));
        otherwise,     n_poles=sum(arrayfun(@(r) numel(r.A),Rfit(:)));
    end
end
fprintf('%s, pole sharing "%s": %d distinct poles in total, fit error %.1f dB.\n',fit_method,pole_sharing,n_poles,Rfit_errordB);

resp=freqresp(Rfit,S.Frequencies);
% Fit error per S-parameter entry in dB (see error_metric). The fit's own
% Rfit_errordB is one value for all entries / per column with shared poles,
% so it cannot rank individual couplings. All-zero entries get -Inf.
Sdata=S.Parameters;
err_num=sum(abs(resp-Sdata).^2,3);
err_den=sum(abs(Sdata).^2,3);
switch error_metric
    case "absolute", Rfit_errordB_mat=10*log10(err_num/numel(S.Frequencies));
    case "relative", Rfit_errordB_mat=10*log10(err_num./err_den);
    otherwise,       error('error_metric must be "absolute" or "relative".');
end
Rfit_errordB_mat(err_den==0)=-Inf;

port_names=get_port_names_from_snp_comments(fullfile(filepath,filename));

%% Plot data vs fit
% S_{a,m} for every antenna port a and every port m to plot. Antenna ports
% are those whose name is not a spherical mode name (TE_L*/TM_L*), e.g.
% "Feed", "load1". One figure per port m, one tile per antenna port;
% magnitude on the left axis, phase on the right; data solid, fit dashed.

is_mode_port=~cellfun(@isempty,regexp(port_names,'^T[EM]_L','once'));
antenna_ports=find(~is_mode_port);
bad=mode_ports_to_plot(mode_ports_to_plot<1 | mode_ports_to_plot>S.NumPorts);
if ~isempty(bad)
    error('mode_ports_to_plot contains port(s) outside 1..%d: %s',S.NumPorts,mat2str(bad));
end
not_mode=mode_ports_to_plot(~is_mode_port(mode_ports_to_plot));
if ~isempty(not_mode)
    warning('mode_ports_to_plot includes non-mode port(s): %s',mat2str(not_mode));
end
antenna_names=string(port_names(antenna_ports));
fprintf('Antenna ports: %s\n',strjoin(compose('%d (%s)',antenna_ports(:),antenna_names(:)),', '));

% Rank every S_{a,m} (a = antenna port, m = any port; all-zero entries
% ignored) by fit error and add the ports m of the n_worst worst ones.
err_am=Rfit_errordB_mat(antenna_ports,:);
[~,order]=sort(err_am(:),'descend');
worst_ports=[];
for k=order(:)'
    [ia,im]=ind2sub(size(err_am),k);
    if isinf(err_am(ia,im)), break, end
    fprintf('Worst-fit coupling #%d: S_{%d,%d} (%s, %s), %.1f dB\n',find(order==k),antenna_ports(ia),im,port_names{antenna_ports(ia)},port_names{im},err_am(ia,im));
    if ~ismember(im,worst_ports), worst_ports(end+1)=im; end
    if numel(worst_ports)>=n_worst, break, end
end
worst_m=worst_ports(1);
plot_ports_always=[];
if always_plot_antenna_ports, plot_ports_always=antenna_ports(:)'; end
plot_mode_ports=unique([mode_ports_to_plot(:)' worst_ports plot_ports_always],'stable');
% Figures ordered from worst to best fit (by the worst coupling of each port).
port_worst_err=arrayfun(@(m) max(Rfit_errordB_mat(antenna_ports,m)),plot_mode_ports);
[~,ord]=sort(port_worst_err,'descend');
plot_mode_ports=plot_mode_ports(ord);
port_worst_err=port_worst_err(ord);
if isempty(mode_ports_to_plot)
    fprintf('mode_ports_to_plot is empty: plotting only the worst-fit port(s) %s.\n',mat2str(plot_mode_ports));
else
    extra=setdiff(plot_mode_ports,mode_ports_to_plot);
    if ~isempty(extra)
        fprintf('Adding port(s) %s to the plots (not in mode_ports_to_plot).\n',mat2str(extra));
    end
end

for pm=1:numel(plot_mode_ports)
    m=plot_mode_ports(pm);
    port_m_escaped=strrep(port_names{m},'_','\_');
    fig_name=sprintf('#%d worst fit: port %d %s (%.1f dB)',pm,m,port_names{m},port_worst_err(pm));
    if ~isempty(mode_ports_to_plot) && ~ismember(m,mode_ports_to_plot), fig_name=fig_name+" (added)"; end
    figure('Name',fig_name,'NumberTitle','off')
    tiledlayout(1,numel(antenna_ports),'TileSpacing','compact','Padding','compact')
    [~,a_ord]=sort(Rfit_errordB_mat(antenna_ports,m),'descend');   % tiles: worst coupling first
    for a=reshape(antenna_ports(a_ord),1,[])
        nexttile
        data_am=rfparam(S,a,m);
        fit_am=squeeze(resp(a,m,:));
        yyaxis left
        h=plot(S.Frequencies,abs(data_am),'-',S.Frequencies,abs(fit_am),'--');
        h(1).LineWidth=1; h(2).LineWidth=1.5;   % thick dashed fit stays visible over the data
        ylabel('Magnitude')
        yyaxis right
        h=plot(S.Frequencies,(180/pi)*angle(data_am),'-',S.Frequencies,(180/pi)*angle(fit_am),'--');
        h(1).LineWidth=1; h(2).LineWidth=2;
        ylabel('Phase (degrees)')
        xlabel('Frequency (Hz)')
        port_a_escaped=strrep(port_names{a},'_','\_');
        title({sprintf('$$S_{%d,%d}$$: %s, %s',a,m,port_a_escaped,port_m_escaped),sprintf('Fit error (%s MSE): %.1f dB',error_metric,Rfit_errordB_mat(a,m))},'Interpreter','latex')
    end
end

%% SPICE netlist and fit export
% "Column" and "Individual" netlists get the pole sharing appended to their name.

ckt_tag="";
if pole_sharing~="All", ckt_tag="_"+pole_sharing; end
ckt_path=fullfile(rfit_dir,filename+ckt_tag+".ckt");
generateSPICE(Rfit,ckt_path,S.Impedance,"SeparateReferenceNodes",true)

% Save the fit next to its netlist for Simulink (matlab/workflows/ook/run_ook_simulink_rf.m,
% rf_fit = "file"). Rfit is a rational object ("rational") or an N x N
% rfmodel.rational array ("rationalfit"); Z0 is the reference impedance.
Z0=S.Impedance;
fit_mat=fullfile(rfit_dir,filename+ckt_tag+"_Rfit.mat");
save(fit_mat,"Rfit","Rfit_errordB","fit_method","pole_sharing","tolerance_dB","max_poles","Z0");
fprintf('Saved the fit (%s, %d poles) to %s\n',fit_method,n_poles,fit_mat);

% Rename the generic generateSPICE node names to the SnP port names (from
% its "! Port[N] Name" comments). Overwrites ckt_path in place.
fix_rational_fit_netlist(ckt_path, fullfile(filepath,filename))

% Add a 1 TOhm gmin resistor from every internal node to ground to help
% convergence (port nodes skipped). Also overwrites ckt_path in place.

add_gmin_to_rational_fit_netlist(ckt_path)
%% Local functions

function [worst_dB,worst_ij]=local_worst_rel_err(R,S)
% Worst relative fit error over the entries, 10*log10(sum|fit-data|^2/sum|data|^2),
% and its [row col] index; all-zero entries are ignored.
H=freqresp(R,S.Frequencies);  D=S.Parameters;
den=sum(abs(D).^2,3);  E=10*log10(sum(abs(H-D).^2,3)./den);  E(den==0)=-Inf;
[worst_dB,k]=max(E(:));  [i,j]=ind2sub(size(E),k);  worst_ij=[i j];
end