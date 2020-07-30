function rxsignal_gen(obj, varargin)

for rx = 1: obj.num_ant
    rxsignal = 0;
    for tx = 1: obj.num_ant
        channel = reshape(obj.time_channel(rx, tx, :), 1, numel(obj.time_channel(rx, tx, :)));
        txsignal = reshape(obj.tx_waveform(tx, :), 1, numel(obj.tx_waveform(tx, :)));
        rxsignal = rxsignal + conv(channel, txsignal);
    end
    obj.rx_waveform(rx, :) = rxsignal;
end



if nargin>1
    SNRdB = varargin{1};
    SNR = 10^(SNRdB/10);   % Linear SNR
    % Calculate noise gain
    N0 = 1/(sqrt(2*obj.num_ant*obj.NFFT)*SNR);
    
    % Create additive white Gaussian noise
    noise = N0*complex(randn(size(obj.rx_waveform)),randn(size(obj.rx_waveform)));
    
    obj.rx_waveform = obj.rx_waveform + noise;
    dbg = 1;
end
dbg = 1;
end

