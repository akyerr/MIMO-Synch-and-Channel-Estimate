classdef OFDM < handle
    % OFDM parameters
    
    properties
        fs = 960e3;
        NFFT = 64;
        CP = 16;
        num_databins = 62;
        num_synchbins = 62;
        bin_spacing = 15e3;
        mod_type = 'QPSK';
        samp_per_symb = [];
        num_ant = 2;
        data_bin_ind = [];
        synch_bin_ind = [];
        num_datasymb = [];
        num_synchsymb = [];
        total_numsymb = [];
    end
    
    methods
        function obj = OFDM(varargin)
            % NFFT, bin_spacing, CP, n_used_databins, mod_type
            for i = 1: nargin
                switch i
                    case 1
                        obj.NFFT = varargin{1};
                    case 2 
                        obj.bin_spacing = varargin{2};
                    case 3
                        obj.CP = varargin{3};
                    case 4
                        obj.num_databins = varargin{4};
                    case 5
                        obj.mod_type = varargin{5};
                    case 6
                        obj.num_ant = varargin{6};
                end
            end
            
            obj.samp_per_symb = obj.NFFT + obj.CP;
            
            obj.data_bin_ind = rem(obj.NFFT+[-obj.num_databins/2:-1,1: obj.num_databins/2], obj.NFFT)+1;
        end
        
        
    end
end

