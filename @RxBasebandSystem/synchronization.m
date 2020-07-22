function synchronization(obj)


% Currently no synchronization - just extract and remove synch symbols
n_total_synch = length(obj.symb_pattern) - sum(obj.symb_pattern);


if n_total_synch == 0
    obj.buffer_data_rx_time = obj.buffer_rx_time;
end


% start position of all synchs
synch_start_pos = 1: sum(obj.synch_data)*obj.samp_per_symb: length(obj.symb_pattern)*obj.samp_per_symb;
for i = 1: n_total_synch
    
    for ant = 1: obj.Nr
        synch_start = synch_start_pos(i);
        synch_fin = synch_start + obj.n_synch*obj.samp_per_symb - 1;
        
        data_start = synch_fin + 1;
        data_fin = data_start + obj.synch_data(2)*obj.samp_per_symb - 1;
        
        start = (i-1)*obj.synch_data(2)*obj.samp_per_symb + 1;
        fin = i*obj.synch_data(2)*obj.samp_per_symb;
        obj.buffer_data_rx_time(ant, start: fin) = obj.buffer_rx_time(ant, data_start: data_fin);
    end
end
dbg = 1;
end

