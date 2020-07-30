function channel_gen(obj)




dbg = 1;

obj.time_channel = zeros(obj.num_ant, obj.num_ant, obj.max_impulse);
obj.freq_channel = zeros(obj.num_ant, obj.num_ant, obj.NFFT);
obj.freq_chan_usedbins = zeros(obj.num_ant, obj.num_ant, obj.num_databins);
h = cell(size(obj.h0));
for rx = 1: obj.num_ant
    for tx = 1: obj.num_ant
        h{rx, tx} = obj.h0{rx, tx};
%         [~, obj.max_tap] = max(h{1, 1});
        obj.time_channel(rx, tx, 1: length(h{rx,tx})) = h{rx,tx}/norm(h{rx,tx});
        obj.freq_channel(rx, tx,:) = fft(obj.time_channel(rx, tx, 1: length(h{rx,tx})), obj.NFFT);
        obj.freq_chan_usedbins(rx,tx,:) = obj.freq_channel(rx, tx, obj.data_bin_ind);
    end
end

dbg = 1;
obj.genie_channel_time = obj.time_channel;
end

