clc
clear
close all
dbstop if error

tic
snr_diag = 100; % dB

scenarios = {'4G5G', 'MIMOTest', 'Fading', 'IndoorA', 960e3*0.9, 15e3, 0;
    '4G5G', 'MIMOTest', 'Fading', 'IndoorB', 960e3*0.9, 15e3, 0};

num_cases = size(scenarios, 1);


for cc = 1: num_cases
    sys_scenario = scenarios{cc, 1};
    MIMO_method = scenarios{cc, 2};
    wireless_channel = scenarios{cc, 3};
    channel_profile = scenarios{cc, 4};
    channel_band = scenarios{cc, 5};
    bin_spacing = scenarios{cc, 6};
    diagnostic = scenarios{cc, 7};
    
    if strcmpi(MIMO_method, 'None')
        num_ant = 1;
    elseif strcmpi(MIMO_method, 'MIMOTest')
        num_ant = 2;
    else
        num_ant = 2; % change this for various MIMO sizes
    end
    
    system = SystemModel(sys_scenario, MIMO_method, num_ant, wireless_channel, channel_profile, ...
        channel_band, bin_spacing, diagnostic, snr_diag);
    
    NFFT = system.NFFT;
    bin_spacing = system.bin_spacing;
    
    if strcmp(system.CP_type, 'Normal')==1
        CP = round(NFFT/4);  %cyclic prefix (IN Samples !!)
    elseif strcmp(system.CP_type,'Extended')==1
        CP = round(NFFT/4 + NFFT/8);  %cyclic prefix (IN Samples !!)
    end
    
    num_databins = 0.75*NFFT;
    num_synchbins = NFFT - 2;
    %     num_databins = 10;
    mod_type = system.mod_type;
    num_ant = system.num_ant;
    bits_per_bin = system.bits_per_bin;
    OFDM_par = OFDM(NFFT, bin_spacing, CP, num_databins, mod_type, num_ant);
    
    
    if diagnostic == 0
        num_iter = length(system.ebno_dB);
    else
        num_iter = 1;
    end
    
    if strcmpi(MIMO_method, 'MIMOTest')
        
        symb_pattern = [1, 0, 2, 0;
                        0, 1, 0, 2]; % 0 - zeros, 1 - synch, 2 - data
        num_symbols = 4; % per subframe? Yes. Per antenna? Yes.
        
        num_data_symb = 1; % per antenna
        
    else
%         symb_pattern here should be automated based on number of subframes and
%         number of symbols per subframe. num_symbols can then simply be
%         the length of symb_pattern.
        error('Currently not supported');
    end
    
    binary_info = randi([0, 1], num_ant, num_databins*bits_per_bin*num_data_symb);
    
    
    Caz = SynchSignal(CP, num_synchbins, num_ant, NFFT);
    ZChu = Caz.ZChu;
    
    OFDM_par.synch_bin_ind = Caz.synch_bin_ind;
    multiant_sys = MultiAntennaSystem(OFDM_par, system, num_symbols);
    
    multiant_sys.multiant_binarymap(symb_pattern, binary_info)
    
end

