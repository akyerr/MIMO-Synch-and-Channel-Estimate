classdef SynchSignal < handle
    properties
        CP = 16;
        num_synchbins = 62;
        num_ant = 2;
        NFFT = 64;
        ZChu = [];
        synch_bin_ind = [];
        time_synchsymb = [];
    end
    methods
        function obj = SynchSignal(varargin)
            for i = 1: nargin
                switch i
                    case 1
                        obj.CP = varargin{1};
                    case 2
                        obj.num_synchbins = varargin{2};
                    case 3
                        obj.num_ant = varargin{3};
                    case 4
                        obj.NFFT = varargin{4};
                end
            end
            
            obj.synch_bin_ind = rem(obj.NFFT+[-obj.num_synchbins/2:-1,1: obj.num_synchbins/2], obj.NFFT)+1;
            
            dbg = 1;
        end
        
        function zadoff_chu_gen(obj, prime_nos)
            
            if length(prime_nos) ~= obj.num_ant
                error('Num of primes has to equal num of antennas')
            end
            
            for ant = 1: obj.num_ant
                
                prime_no = prime_nos(ant);
                seq_len = obj.num_synchbins;
                %Zadoff Chu CAZAC sequence
                switch rem(seq_len, 2)
                    case 0
                        obj.ZChu(ant, :) = exp(-1i*(2*pi/seq_len)*prime_no*((0:seq_len-1).^2)/2);
                    case 1
                        obj.ZChu(ant, :) = exp(-1i*(2*pi/seq_len)*prime_no*((0:seq_len-1).*(1:seq_len))/2);
                end
                
                
                
                freq_synchsymb = zeros(1, obj.NFFT);
                freq_synchsymb(obj.synch_bin_ind) = obj.ZChu(ant, :);
                
                synch_ifftout = ifft(freq_synchsymb, obj.NFFT);
                synch_cp = synch_ifftout(end - obj.CP + 1: end);
                obj.time_synchsymb(ant, :) = [synch_cp, synch_ifftout];
            end
        end
        
    end
end