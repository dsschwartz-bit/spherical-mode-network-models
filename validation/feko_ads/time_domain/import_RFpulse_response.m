%% import_RFpulse_response.m
% Imports the Tx-loop / Rx-bowtie RF-pulse transient exported from ADS (Tx_loop_Rx_bowtie_RFpulse.csv).
%
% Input:  csvFile - path of the CSV, set by the calling script
%         (default: Tx_loop_Rx_bowtie_RFpulse.csv next to this script).
% Output (workspace variables, same names as import_impulse_response.m):
%   time                  [s]  (transient time steps, not uniform)
%   ADS_RationalFitModes  [V]  node FeedRx   - Rx feed of the rational-fit modal blocks
%   FEKO_S2P              [V]  node FeedRx_2 - Rx feed of the FEKO .s2p SnP block
%   Vgen                  [V]  node Vgen     - generator voltage
%
% Column mapping: in schematic Txloop_Rxbowtie_internal_ports_transient_ARRANGED, net FeedRx_2 is
% at the FEKO .s2p SnP block and FeedRx at the rational-fit bowtie (s_dom_q_bowtie_..._s6p).
% readADSTimeCsv checks the header references, so an export with a different column order stops
% with an error.

if ~exist('csvFile', 'var')
    csvFile = fullfile(fileparts(mfilename('fullpath')), "Tx_loop_Rx_bowtie_RFpulse.csv");
end
expectedRefs = ["FeedRx", "FeedRx_2", "Vgen"];
varNames     = ["time", "ADS_RationalFitModes", "FEKO_S2P", "Vgen"];   % file column order

[time, cols] = readADSTimeCsv(csvFile, expectedRefs);
ADS_RationalFitModes = cols(:, 1);
FEKO_S2P             = cols(:, 2);
Vgen                 = cols(:, 3);
clear cols csvFile expectedRefs varNames
