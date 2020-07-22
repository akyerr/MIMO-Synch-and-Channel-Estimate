classdef SystemModel < handle
    properties
        sys_scenario = 'WIFI';
        diagnostic = 0;
        wireless_channel = 'Fading';
        channel_profile = 'IndoorA';
        channel_band = 20e6;
        bin_spacing = 15e3;
        case_number = 1;
        signal_dtype = 'Complex';
        NFFT = [];
        num_synch_bins = [];
        mod_type = 'QPSK';
        bits_per_bin = 2;
        PHYchan = 'Data';
        fs = [];
        CP_type = 'Normal';
        num_ant = 2;
        param_est = 'Ideal';
        MIMO_method = 'SpMux';
        stream_size = 1;
        ebno_dB = [];
        num_subframes = [];
        synch_data = [1, 2]; % 1 synchs for every 2 data symbols
    end
    methods
        function obj = SystemModel(varargin)
            for i = 1: nargin
                switch i
                    case 1
                        obj.sys_scenario = varargin{1};
                    case 2
                        obj.MIMO_method = varargin{2};
                    case 3
                        obj.num_ant = varargin{3};
                    case 4
                        obj.wireless_channel = varargin{4};
                    case 5
                        obj.channel_profile = varargin{5};
                    case 6
                        obj.channel_band = varargin{6};
                    case 7
                        obj.bin_spacing = varargin{7};
                    case 8
                        obj.diagnostic = varargin{8};
                end
            end
            
            
            obj.NFFT = 2^(ceil(log2(round(obj.channel_band/obj.bin_spacing))));  %4G/5G
            obj.num_synch_bins = obj.NFFT-2;
            
            obj.signal_dtype = 'Complex'; %Most of the time assume complex, affects AWGN addition
            obj.PHYchan = 'Data'; %'Data', 'Synch', 'Security':  Security is not used
            obj.mod_type = 'QPSK'; % per bin
            obj.bits_per_bin = 2;
            obj.fs = obj.bin_spacing*obj.NFFT;
            
            if obj.diagnostic == 0
                obj.ebno_dB = [6, 7, 8, 9, 10, 14, 16, 20, 24];
            else
                obj.ebno_dB = snr_diag;
            end
            
            obj.CP_type = 'Normal'; % 'Normal' or 'Extended'
            obj.param_est = 'Ideal'; % 'Ideal' or 'Estimated
            obj.stream_size = obj.num_ant;
        end
    end
end