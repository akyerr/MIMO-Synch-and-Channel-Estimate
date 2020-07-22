classdef MultiAntennaSystem < handle
    
    properties
        num_ant = 2;
        OFDM_par = [];
        system = [];
        NFFT = [];
        CP = [];
        max_impulse = [];
        time_channel = [];
        freq_channel = [];
        freq_chan_usedbins = [];
        data_bin_ind = [];
        num_databins = [];
        genie_channel_time = [];
        max_tap = [];
        epsval = 1e-5;
        fs = [];
        channel_profile = [];
        h0 = [];
        seed0 = [];
        
        tx_symbs = [];
        rx_symbs = [];
        
        tx_waveform = [];
        rx_waveform = [];
        num_symbols = [];
        samp_per_symb = [];
    end
    
    methods
        function obj = MultiAntennaSystem(varargin)
            for i = 1: nargin
                switch i
                    case 1
                        obj.OFDM_par = varargin{1};
                        obj.NFFT = obj.OFDM_par.NFFT;
                        obj.CP = obj.OFDM_par.CP;
                        obj.data_bin_ind = obj.OFDM_par.data_bin_ind;
                        obj.num_databins = obj.OFDM_par.num_databins;
                        obj.num_ant = obj.OFDM_par.num_ant;
                        obj.samp_per_symb = obj.OFDM_par.samp_per_symb;
                    case 2
                        obj.system = varargin{2};
                        obj.fs = obj.system.fs;
                        obj.channel_profile = obj.system.channel_profile;
                    case 3
                        obj.num_symbols = varargin{3};
                end
            end
            
            
            obj.max_impulse = obj.NFFT;
            
            obj.tx_symbs = zeros(obj.num_ant, obj.NFFT*obj.num_symbols);
            obj.rx_symbs = zeros(obj.num_ant, obj.NFFT*obj.num_symbols);
            
            obj.tx_waveform = zeros(obj.num_ant, obj.samp_per_symb*obj.num_symbols);
            obj.tx_waveform = zeros(obj.num_ant, obj.samp_per_symb*obj.num_symbols + obj.max_impulse - 1);
            
            
            %% channels
            
            [PDPdB, PathDelay] = ChannelProfile(obj.channel_profile);
            test_case = 0; % 0 - LTE channel, 1 - test channel 1
            if obj.num_ant == 1
                if test_case == 0
                    
                    rayChan = comm.RayleighChannel(...
                        'SampleRate', obj.fs, ...
                        'PathDelays', PathDelay, ...
                        'AveragePathGains',PDPdB, ...
                        'NormalizePathGains',true, ...
                        'MaximumDopplerShift',30, ...
                        'DopplerSpectrum',{doppler('Gaussian',0.6)}, ...
                        'RandomStream','mt19937ar with seed', ...
                        'Seed',obj.seed0, ...
                        'PathGainsOutputPort',true,...
                        'Visualization','off');
                    hprime=rayChan([1;zeros(obj.NFFT,1)]); %MIMO Columns
                    
                    NW = 5;
                    L=1;
                    while L < length(hprime) && sum(abs(hprime(L:L+NW-1).*conj(hprime(L:L+NW-1))))/NW > obj.epsval
                        L = L+1;
                    end
                    backI = 0;
                    while hprime(L-backI) == 0
                        backI = backI+1;
                    end
                    L = L-backI;
                    
                    obj.h0{1,1}=hprime(1:L);
                    
                    %                     obj.h0{1, 1}
                    %                     obj.max_tap
                    dbg = 1;
                elseif test_case == 1
                    obj.h0{1,1}=[0.3977,  -0.1988,  0.0994,  -0.0398, 0.7954 - 0.3977i].';
                    [~, obj.max_tap] = max(obj.h0{1, 1});
                    
                end
            elseif obj.num_ant == 2
                if test_case == 1
                    obj.h0{1,1}=[0.3977,  0.7954 - 0.3977i,  -0.1988,  0.0994,  -0.0398].';
                    obj.h0{1,2}=[0.8423i,  0.5391, 0, 0, 0].';
                    obj.h0{2,1}=[0.1631,  -0.0815 + 0.9784i,  0.0978, 0, 0].';
                    obj.h0{2,2}=[0.0572i,  0.3659i, 0.5717 - 0.5717i, 0.4574, 0].';
                end
            else
                error('currently not supported')
            end
            
            
        end
        
        
    end
end

