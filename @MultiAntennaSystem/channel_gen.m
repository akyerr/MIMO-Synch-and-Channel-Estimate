function channel_gen(obj)

obj.time_channel = zeros(obj.num_ant, obj.num_ant, obj.max_impulse);
obj.freq_channel = zeros(obj.num_ant, obj.num_ant, obj.NFFT);
obj.freq_chan_usedbins = zeros(obj.num_ant, obj.num_ant, obj.num_synchbins);
% obj.freq_chan_usedbins_synch = zeros(obj.num_ant, obj.num_ant, obj.num_synchbins);
h = cell(size(obj.h0));
chan_seenby_rx = cell(obj.num_ant, obj.num_ant);
for rx = 1: obj.num_ant
    for tx = 1: obj.num_ant
        h{rx, tx} = obj.h0{rx, tx};
        obj.time_channel(rx, tx, 1: length(h{rx,tx})) = h{rx,tx}/norm(h{rx,tx});
        obj.freq_channel(rx, tx,:) = fft(obj.time_channel(rx, tx, 1: length(h{rx,tx})), obj.NFFT);
        obj.freq_chan_usedbins(rx, tx, :) = obj.freq_channel(rx, tx, obj.synch_bin_ind);
        
%         chan_seenby_rx0 = reshape(obj.time_channel(rx, tx, :), 1, numel(obj.time_channel(rx, tx, :)));
%         ind = 10000;
%         
%         while ind ~= 1
%             chan_seenby_rx0 = circshift(chan_seenby_rx0, 1);
%             [~, ind] = max(abs(chan_seenby_rx0));
%             dbg = 1;
%         end
%         chan_seenby_rx{rx, tx} = chan_seenby_rx0;
%         dbg = 1;
    end
end

dbg = 1;
obj.genie_channel_time = obj.time_channel;
end

