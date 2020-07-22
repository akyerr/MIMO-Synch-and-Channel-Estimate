function tx_waveform = OFDM_mod(obj, tx_vector, ZChu, synch_data)
% number of data symbols in current subframe
obj.num_datasymb = size(tx_vector, 2)/obj.num_databins;

% 0 - synch, 1 - data
symb_pattern0 = [zeros(1, synch_data(1)), ones(1, synch_data(2))];
symb_pattern = repmat(symb_pattern0, 1, ceil(obj.num_datasymb/synch_data(2)));

if sum(symb_pattern) > obj.num_datasymb
    symb_pattern = symb_pattern(1: end - (sum(symb_pattern)-obj.num_datasymb));
end

obj.num_synchsymb = length(find(symb_pattern==0));
obj.total_numsymb = length(symb_pattern);




data_time = zeros(obj.num_ant, obj.num_datasymb*obj.samp_per_symb);
for ant = 1: obj.num_ant
    for symb = 1: obj.num_datasymb
        freq_data = zeros(1, obj.NFFT);
        fbin_start = (symb-1)*obj.num_databins + 1;
        fbin_end = fbin_start + (obj.num_databins -1);
        
        freq_data(obj.data_bin_ind) = tx_vector(ant, fbin_start: fbin_end);
        
        data_ifft = ifft(freq_data, obj.NFFT);
        data_cp = data_ifft(end-obj.CP+1: end); % CP
        
        time_start = (symb-1)*obj.samp_per_symb + 1;
        time_end = time_start + (obj.samp_per_symb-1);
        
        data_time(ant, time_start: time_end) = [data_cp, data_ifft]; % Add CP
        
    end
end

freq_synch = zeros(1, obj.NFFT);
freq_synch(obj.synch_bin_ind) = ZChu;

synch_ifft = ifft(freq_synch);
synch_cp = synch_ifft(end-obj.CP+1: end);
synch_time0 = [synch_cp, synch_ifft];

synch_time = zeros(obj.num_ant, obj.num_synchsymb*obj.samp_per_symb);
synch_time(1, :) = repmat(synch_time0, 1, obj.num_synchsymb);

tx_waveform = zeros(obj.num_ant, obj.total_numsymb*obj.samp_per_symb);

total_symbcount = 0;
synch_symbcount = 0;
data_symbcount = 0;
for symb = symb_pattern
    symb_start = total_symbcount*obj.samp_per_symb + 1;
    symb_end = symb_start + (obj.samp_per_symb-1);
    if symb == 0
        synch_start = synch_symbcount*obj.samp_per_symb + 1;
        synch_end = synch_start + (obj.samp_per_symb-1);
        
        tx_waveform(:, symb_start: symb_end) = synch_time(:, synch_start: synch_end);
        synch_symbcount = synch_symbcount + 1;
    elseif symb == 1
        data_start = data_symbcount*obj.samp_per_symb + 1;
        data_end = data_start + (obj.samp_per_symb-1);
        
        tx_waveform(:, symb_start: symb_end) = data_time(:, data_start: data_end);
        data_symbcount = data_symbcount + 1;  
    end
    
    total_symbcount = total_symbcount + 1;
end

dbg = 1;

end

