function synchronize(obj, synch_data_pattern, symb_pattern)
CP = obj.CP;
NFFT = obj.NFFT;

num_tx_symb = size(symb_pattern, 2);
max_exp_data_elements = num_tx_symb*obj.samp_per_symb;


start_samp = 4;

buffer_rx_time = obj.rx_waveform(:, start_samp: end);
mask = obj.synch_ref_freq;

rx_power = power_estimate(buffer_rx_time(1, :)); % power check only on antenna 1

if rx_power > obj.power_requirements
    
    est_chan_avg = cell(obj.num_ant, obj.num_ant);
    
    del = obj.max_tap;
    for rx_ant = 1: obj.num_ant
        for tx_ant = 1: obj.num_ant
            
            disp(['At rx ant ', num2str(rx_ant), ' from tx ant ', num2str(tx_ant)]);
            
            [max_corr_val, max_corr_ind] = correlate(buffer_rx_time(rx_ant, :), mask(tx_ant, :), NFFT, CP);
            
            obs_synch_freq = fft(buffer_rx_time(rx_ant, max_corr_ind: max_corr_ind + NFFT), NFFT);
            
            
            obs_synch_at_usedbins = obs_synch_freq(obj.synch_bin_ind);

            synch_ref = obj.ZChu(tx_ant, :);
            
            est_chan_freq = (obs_synch_at_usedbins.*conj(synch_ref))./abs(synch_ref);
            
            est_chan_pow = sum(est_chan_freq.*conj(est_chan_freq))/numel(est_chan_freq);
            
            est_chan_freq = est_chan_freq*sqrt(1/(est_chan_pow));
            actual_chan_freq = reshape(obj.freq_chan_usedbins(rx_ant, tx_ant, :), 1, numel(obj.freq_chan_usedbins(rx_ant, tx_ant, :)));
            
            est_chan_allbins = zeros(1, NFFT);
            actual_chan_freq_allbins = zeros(1, NFFT);
            est_chan_allbins(obj.synch_bin_ind) = est_chan_freq;
            actual_chan_freq_allbins(obj.synch_bin_ind) = actual_chan_freq;
            
            
            chan_seenby_rx1 = reshape(obj.genie_channel_time(rx_ant, tx_ant, :), 1, numel(obj.genie_channel_time(rx_ant, tx_ant, :)));
%             ind = 10000;
%             
%             while ind ~= 1
%                 chan_seenby_rx1 = circshift(chan_seenby_rx1, 1);
%                 [~, ind] = max(abs(chan_seenby_rx1));
%                 dbg = 1;
%             end
            
            figure()
            xax = (0:NFFT-1)*obj.fs/obj.NFFT;
            yax1 = 20*log10(abs(ifft(est_chan_allbins, NFFT)));
            yax2 = 20*log10(abs(chan_seenby_rx1));
            %             yax3 = 20*log10(abs(ifft(actual_chan_freq_allbins, NFFT)));
            title('Time Domain')
            plot(xax, yax1, 'r', xax, yax2,'b')
            legend({'Estimated', 'Actual rotated'})
            
%             figure()
%             xax = (0:obj.NFFT-1)*obj.fs/obj.NFFT;
% %             yax1 = 20*log10(abs(est_chan_freq));
% %             yax2 = 20*log10(abs(actual_chan_freq));
%             yax1 = 20*log10(abs(est_chan_allbins));
%             yax2 = 20*log10(abs(actual_chan_freq_allbins));
%             title('Frequency Domain')
%             plot(xax, yax1, 'r', xax, yax2,'b')
%             legend({'Estimated', 'Actual'})
            dbg = 1; %#ok
            
        end
    end
else
    error('rx signal is not strong enough');
end




dbg = 10;
end


function sig_power = power_estimate(input_data)
sig_power = sum(input_data.*conj(input_data))/numel(input_data);
end

function [max_corr_val, max_corr_ind] = correlate(rx_signal, mask, NFFT, CP)

window_size = length(mask);
window_slide_dist = 1;
num_windows = ceil(((length(rx_signal) - window_size)/window_slide_dist) + 1);

corr_val = zeros(1, num_windows);

% autocorr in Taylor python code
for win = 1: num_windows
    
    start = (win-1)*window_slide_dist +1;
    fin = start + window_size - 1;
    
    rx_signal_freq = fft(rx_signal(start: fin), NFFT);
    
    corr_val(win) = abs(sum(rx_signal_freq*mask'));
    
end

% figure()
% plot(corr_val)
% xlabel('Window index')
% ylabel('Correlation value')

[max_corr_val, max_corr_ind] = max(corr_val);


dbg = 1;
end






