classdef RxBasebandSystem < handle
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        num_ant = 2;
        NFFT = [];
        CP = [];
        samp_per_symb = [];
        data_bin_ind = [];
        synch_bin_ind = [];
        freq_channel = [];
        time_channel = [];
        genie_channel_time = [];
        multiant_sys = [];
        Caz = [];
        synch_ref_freq = [];
        rx_waveform = [];
        system = [];
        stride_val = [];
        time_synch_ref = [];
        corr_obs = [];
        start_samp = [];
        del_mat = [];
        diagnostic = [];
        fs = [];
        OFDM_par = [];
        num_synchsymb = [];
        est_chan_freqP = [];
        est_chan_freqN = [];
        est_chan_time = [];
        est_synch_freq = [];
        samp_delay = [];
        max_tap = [];
        synch_ref_time = [];
        tx_waveform = [];
        power_requirements = 0;
        freq_chan_usedbins = [];
        ZChu = [];
    end
    
    methods
        function obj = RxBasebandSystem(varargin)
            for i = 1: nargin
                switch i
                    case 1
                        obj.multiant_sys = varargin{1};
                        obj.num_ant = obj.multiant_sys.num_ant;
                        obj.NFFT = obj.multiant_sys.NFFT;
                        obj.CP = obj.multiant_sys.CP;
                        obj.samp_per_symb = obj.multiant_sys.OFDM_par.samp_per_symb;
                        obj.data_bin_ind = obj.multiant_sys.data_bin_ind;
                        obj.freq_channel = obj.multiant_sys.freq_channel;
                        obj.genie_channel_time = obj.multiant_sys.genie_channel_time;
                        
                        obj.max_tap = obj.multiant_sys.max_tap;
                        obj.rx_waveform = obj.multiant_sys.rx_waveform;
                        obj.tx_waveform = obj.multiant_sys.tx_waveform;
                        obj.freq_chan_usedbins = obj.multiant_sys.freq_chan_usedbins;
                    case 2
                        obj.Caz = varargin{2};
                        obj.synch_bin_ind = obj.Caz.synch_bin_ind;
                        obj.synch_ref_freq = obj.Caz.freq_synchsymb;
                        obj.ZChu = obj.Caz.ZChu;
                        obj.synch_ref_time = obj.Caz.time_synchsymb;
                    case 3
                        obj.system = varargin{3};
                        obj.diagnostic = obj.system.diagnostic;
                        obj.fs = obj.system.fs;
                    case 4
                        obj.OFDM_par = varargin{4};
                end
            end
        end
        

    end
end

