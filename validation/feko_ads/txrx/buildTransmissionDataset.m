function d = buildTransmissionDataset(filepath, fekoMagDownsampleHz, fekoPhaseDownsampleHz)
%BUILDTRANSMISSIONDATASET  Build one ADS/FEKO transmission dataset from an .s4p file.
%   d = BUILDTRANSMISSIONDATASET(filepath, fekoMagDownsampleHz, fekoPhaseDownsampleHz)
%
%   Reads the Touchstone file with sparameters/rfparam (RF Toolbox). Ports
%   1-2 are the ADS Tx/Rx pair and ports 3-4 the FEKO Tx/Rx pair (two
%   independent links in one 4-port file; S43/S34/S33/S44 equal the FEKO
%   .s2p S21/S12/S11/S22), so transmission is S21 (ADS) and S43 (FEKO).
%
%   The FEKO magnitude and phase are thinned independently (downsampleFeko)
%   with their own frequency steps [Hz].
%
%   Output struct fields:
%     theta               - theta angle string, parsed from the
%                            filename (the "th<digits>" token)
%     phi                 - phi angle string, parsed from the
%                            filename (the "ph<digits>" token)
%     freq                - dense frequency vector [Hz] (full resolution)
%     ads.mag             - dense ADS |S21| in dB
%     ads.phase           - dense ADS phase(S21) in degrees
%     fekoDS.mag.freq     - FEKO frequency vector [Hz], downsampled to
%                            fekoMagDownsampleHz
%     fekoDS.mag.data     - downsampled FEKO |S43| in dB
%     fekoDS.phase.freq   - FEKO frequency vector [Hz], downsampled to
%                            fekoPhaseDownsampleHz
%     fekoDS.phase.data   - downsampled FEKO phase(S43) in degrees

    s_obj = sparameters(filepath);
    freq = s_obj.Frequencies;

    S21 = rfparam(s_obj, 2, 1);   % ADS Rx <- ADS Tx
    S43 = rfparam(s_obj, 4, 3);   % FEKO Rx <- FEKO Tx

    fekoMagDb    = 20*log10(abs(S43));
    fekoPhaseDeg = angle(S43) * 180/pi;
    adsMagDb    = 20*log10(abs(S21));
    adsPhaseDeg = angle(S21) * 180/pi;

    [~, name, ~] = fileparts(filepath);
    thetaTok = regexp(name, "th(\d+)", "tokens", "once");
    theta = thetaTok{1};
    phiTok = regexp(name, "ph(\d+)", "tokens", "once");
    phi = phiTok{1};

    fekoMagStruct   = downsampleFeko(struct("freq",freq, "data",fekoMagDb,   "tags","x"), fekoMagDownsampleHz);
    fekoPhaseStruct = downsampleFeko(struct("freq",freq, "data",fekoPhaseDeg,"tags","x"), fekoPhaseDownsampleHz);

    d.theta = theta;
    d.phi   = phi;
    d.freq  = freq;
    d.ads.mag   = adsMagDb;
    d.ads.phase = adsPhaseDeg;
    d.fekoDS.mag.freq    = fekoMagStruct.freq;
    d.fekoDS.mag.data    = fekoMagStruct.data;
    d.fekoDS.phase.freq  = fekoPhaseStruct.freq;
    d.fekoDS.phase.data  = fekoPhaseStruct.data;
end