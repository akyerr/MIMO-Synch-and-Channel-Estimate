function rx_waveform = rxsignal_gen(obj, tx_waveform, varargin)

rx_waveform = zeros(obj.num_ant, size(tx_waveform, 2) + obj.max_impulse - 1);
for rx = 1: obj.num_ant
    rxsignal = 0;
    for tx = 1: obj.num_ant
        channel = reshape(obj.time_channel(rx, tx, :), 1, numel(obj.time_channel(rx, tx, :)));
        txsignal = reshape(tx_waveform(tx, :), 1, numel(tx_waveform(tx, :)));
        rxsignal = rxsignal + conv(channel, txsignal);
    end
    rx_waveform(rx, :) = rxsignal;
end



if nargin>2
    SNRdB = varargin{1};
    SNR = 10^(SNRdB/10);   % Linear SNR
    % Calculate noise gain
    N0 = 1/(sqrt(2*obj.num_ant*obj.NFFT)*SNR);
    
    % Create additive white Gaussian noise
    noise = N0*complex(randn(size(rx_waveform)),randn(size(rx_waveform)));
    
    rx_waveform = rx_waveform + noise;
    dbg = 1;
end
dbg = 1;
end

