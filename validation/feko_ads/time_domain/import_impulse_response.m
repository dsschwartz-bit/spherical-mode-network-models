%% import_impulse_response.m
% Imports the Tx-loop / Rx-bowtie impulse responses exported from ADS (Tx_loop_Rx_bowtie_ts.csv).
%
% Input:  csvFile - path of the CSV, set by the calling script
%         (default: Tx_loop_Rx_bowtie_ts.csv next to this script).
% Output (workspace variables, same names as import_RFpulse_response.m):
%   time                  [s]
%   ADS_SnPModes          ADS modal SnP (Touchstone) blocks  ts(Txloop_Rxbowtie_internal_ports_transient..S(2,1))
%   FEKO_S2P              FEKO .s2p placed as an SnP block   ts(S(4,3))
%   ADS_RationalFitModes  ADS rational-fit modal blocks      ts(Txloop_Rxbowtie_internal_ports_transient_ARRANGED..S(2,1))
%
% Column mapping: schematic Txloop_Rxbowtie_internal_ports_transient uses the modal SnP
% (s_dom_q_*_TOUCHESTONE) subcircuits; the _ARRANGED schematic uses the rational-fit .net models
% (s_dom_q_*.net) and holds the FEKO block (Tx_loop_Rx_bowtie_ph90th60_SParameter.s2p).
% readADSTimeCsv checks the header references, so an export with a different column order stops
% with an error.

if ~exist('csvFile', 'var')
    csvFile = fullfile(fileparts(mfilename('fullpath')), "Tx_loop_Rx_bowtie_ts.csv");
end
expectedRefs = ["ts(Txloop_Rxbowtie_internal_ports_transient..S(2,1))", "ts(S(4,3))", ...
                "ts(Txloop_Rxbowtie_internal_ports_transient_ARRANGED..S(2,1))"];
varNames     = ["time", "ADS_SnPModes", "FEKO_S2P", "ADS_RationalFitModes"];   % file column order

[time, cols] = readADSTimeCsv(csvFile, expectedRefs);
ADS_SnPModes         = cols(:, 1);
FEKO_S2P             = cols(:, 2);
ADS_RationalFitModes = cols(:, 3);
clear cols csvFile expectedRefs varNames
