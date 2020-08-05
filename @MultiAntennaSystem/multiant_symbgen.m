function multiant_symbgen(obj, num_symbols)

% num_symbols - number of symbols per antenna

min_pow = 1e-30;

for symb = 1: num_symbols
    scale_factor = 1;
    pow = 0;
    
    for ant = 1: obj.num_ant
        ifft_in = obj.tx_symbs(ant, (symb-1)*obj.NFFT+1: symb*obj.NFFT);
        ifft_out = ifft(ifft_in, obj.NFFT);
        cyclic_prefix = ifft_out(end - obj.CP + 1: end);
        time_symb = [cyclic_prefix, ifft_out];
%         if sum(time_symb*time_symb') > min_pow && ant == 1
        if sum(time_symb*time_symb') > min_pow 
            scale_factor = sqrt(length(time_symb)/sum(time_symb*time_symb'));
        end
        time_symb = time_symb*scale_factor;
        
        obj.tx_waveform(ant, (symb-1)*obj.samp_per_symb+1: symb*obj.samp_per_symb) = time_symb;
        
        pow = pow + var(time_symb);
        
    end
    
    for ant = 1: obj.num_ant
        obj.tx_waveform(ant, (symb-1)*obj.samp_per_symb+1: symb*obj.samp_per_symb) = ...
            obj.tx_waveform(ant, (symb-1)*obj.samp_per_symb+1: symb*obj.samp_per_symb)*(1/sqrt(pow));
    end
end
end