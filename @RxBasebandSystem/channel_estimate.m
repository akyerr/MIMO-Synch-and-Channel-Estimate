function channel_estimate(obj, pilot_perc)

dbg = 1;


for symb = 1: sum(obj.symb_pattern)
    
    
    for ant = 1: obj.Nr
        start = (symb - 1)*obj.samp_per_symb + 1;
        fin = symb*obj.samp_per_symb;
        data_symb = obj.buffer_data_rx_time(start: fin);
        
        data_symb_without_cp = data_symb(obj.CP + 1: end);
        
        data_freq = fft(data_symb_without_cp, obj.NFFT);
        
        scale_factor = length(data_freq)/sum(data_freq.*conj(data_freq));
        
        data_freq = data_freq* sqrt(scale_factor);
        
        f =(-obj.fs/2: obj.delta_f: obj.fs/2 - obj.delta_f);
        figure(10)
        stem(f, abs(data_freq));
        title('Rx data in Frequency Domain - Magnitude')
        
        pilots_freq = data_freq(obj.pilot_bins{ant});
        
        least_sq_est = pilots_freq./obj.pilots;
        
        interpNum = 1/pilot_perc;
        ls_upsamp = upsample(least_sq_est, interpNum);
        
        fs1 = obj.fs*interpNum;
        N1 = obj.NFFT* interpNum;
        indshift1 = [N1/2+1:N1, 1: N1/2];
        s3 = 20 * log10(abs(fft(s2t, N1)));
        xax3 = (-N1/2:N1/2-1)*fs1/N1;
        figure(pltnum)
        pltnum = pltnum +1;
        plot(xax3, s3(indshift1));
        
        
        
        dbg = 1;
        
    end
    dbg = 1;
end

end

