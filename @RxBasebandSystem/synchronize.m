function synchronize(obj, synch_data_pattern, symb_pattern)

num_tx_symb = size(symb_pattern, 2);
max_exp_data_elements = num_tx_symb*obj.samp_per_symb;


start_samp = 4;

buffer_rx_time = obj.rx_waveform(:, start_samp: end);
mask = obj.synch_ref_time;

% figure()
% plot(real(fft(mask(1, obj.CP+1: end))), imag(fft(mask(1, obj.CP+1: end))), '.')
% title('Reference synch FFT')

% Plotting rx waveforms
% for rx_ant = 1: obj.num_ant
%     figure()
%     xax = 1: length(obj.rx_waveform(rx_ant, :));
%     plot(xax, real(obj.rx_waveform(rx_ant, :)), xax, imag(obj.rx_waveform(rx_ant, :)));
%     xlabel('Time')
%     ylabel('Amplitude')
%     title(['Antenna ', num2str(rx_ant), ' Amplitude of Rx Waveform'])
% end


rx_power = power_estimate(buffer_rx_time(1, :)); % power check only on antenna 1

if rx_power > obj.power_requirements
    
    for rx_ant = 1: obj.num_ant
        disp(rx_ant)
        
        for tx_ant = 1: obj.num_ant
            [~, max_corr_ind] = correlate(buffer_rx_time(rx_ant, :), mask(tx_ant, :));
            
            power_est = power_estimate(buffer_rx_time(rx_ant, 1:max_corr_ind));
            
            %% Roadmap of rx symbols and their types
            if max_corr_ind > 1
                if power_est > obj.power_requirements
                    
                    %% Symbols BEFORE synchpt
                    num_full_symbs_before_synchpt = floor((max_corr_ind-1)/(obj.samp_per_symb));
                    num_part_symbs_before_synchpt = ceil((max_corr_ind-1)/(obj.samp_per_symb));
                    
                    txsymb_pattern_before = repmat(synch_data_pattern(tx_ant, :), 1, ceil(num_part_symbs_before_synchpt/length(synch_data_pattern(tx_ant, :))));
                    txsymb_pattern_before = txsymb_pattern_before(1: num_part_symbs_before_synchpt);
                    rxsymb_pattern_before = txsymb_pattern_before(end - num_full_symbs_before_synchpt + 1: end);
                    
                    symb_start_before = zeros(1, num_full_symbs_before_synchpt);
                    symb_type_before = zeros(1, num_full_symbs_before_synchpt);
                    start_ind_before = max_corr_ind;
                    for i = 1: num_full_symbs_before_synchpt
                        start_ind_before = start_ind_before - obj.samp_per_symb;
                        symb_start_before(i) = start_ind_before;
                        symb_type_before(i) = rxsymb_pattern_before(i);
                    end
                    symb_start_before = fliplr(symb_start_before);
                    symb_type_before = symb_type_before; %#ok
                    
                    %% Symbols AFTER synchpt - INCLUDES the synch symbol that was synched to
 
                    num_full_symbs_after_synchpt = floor((max_exp_data_elements - (max_corr_ind-1))/(obj.samp_per_symb)); % Includes the synch that was synched to
                    
                    txsymb_pattern_after = repmat(synch_data_pattern(tx_ant, :), 1, ceil(num_full_symbs_after_synchpt/length(synch_data_pattern(tx_ant, :))));
                    rxsymb_pattern_after = txsymb_pattern_after(end - num_full_symbs_after_synchpt + 1: end);
                    
                    symb_start_after = zeros(1, num_full_symbs_after_synchpt);
                    symb_type_after = zeros(1, num_full_symbs_after_synchpt);
                    start_ind_after = max_corr_ind;
                    
                    for j = 1: num_full_symbs_after_synchpt
                        symb_start_after(j) = start_ind_after;
                        symb_type_after(j) = rxsymb_pattern_after(j);
                        start_ind_after = start_ind_after + obj.samp_per_symb;
                    end
                    
                    
                    rx_symb_start = [symb_start_before, symb_start_after];
                    rx_symb_type = [symb_type_before, symb_type_after];
                    
                    dbg = 1;
                end
            else
                if power_est > obj.power_requirements
                    rx_symb_start = 1: obj.samp_per_symb: max_exp_data_elements;
                    rx_symb_type = symb_pattern(rx_ant, :);
                end
            end
            
            synch_ind = find(rx_symb_type==1);
            if rx_ant == 1
                synch_start = rx_symb_start(synch_ind(1: end));
            else
                synch_start = rx_symb_start(synch_ind(2: end));
            end
            
            synch_end = synch_start + obj.samp_per_symb-1;
            
            %% Channel Estimation
            est_chan_symb = zeros(length(synch_start), obj.NFFT);
            
            for synch = 1: length(synch_start)
                synch_symb = buffer_rx_time(rx_ant, synch_start(synch): synch_end(synch));
                synch_symb_without_cp = synch_symb(obj.CP + 1: end);
                
                obs_synch_freq = fft(synch_symb_without_cp, obj.NFFT);
                
                obs_synch_pow = power_estimate(obs_synch_freq);
                
                obs_synch_freq = (1/sqrt(obs_synch_pow))*obs_synch_freq;
                
                dbg = 1; %#ok
                
                obs_synch_at_usedbins = obs_synch_freq(obj.synch_bin_ind);
                
                synch_ref = obj.synch_ref_freq(tx_ant, :);
                
                est_chan_freq = (obs_synch_at_usedbins.*conj(synch_ref))./abs(synch_ref);
                
                est_chan_pow = sum(est_chan_freq.*conj(est_chan_freq))/numel(est_chan_freq);
                
                est_chan_freq = est_chan_freq*sqrt(1/(est_chan_pow));
                
                
                est_chan_symb(synch, obj.synch_bin_ind) = est_chan_freq;
            end
            
            est_chan_avg = sum(est_chan_symb)/size(est_chan_symb, 1);
            
            genie_channel = reshape(obj.genie_channel_time(rx_ant, tx_ant, :), 1, numel(obj.genie_channel_time(rx_ant, tx_ant, :)));
            figure()
            xax = (0:obj.NFFT-1)*obj.fs/obj.NFFT;
            yax1 = 20*log10(abs(est_chan_avg));
            yax2 = 20*log10(abs(fft(genie_channel, obj.NFFT)));
            plot(xax, yax1, 'r', xax, yax2,'b')
            legend({'Estimated', 'Actual'})
            dbg = 1; %#ok

        end
        dbg = 1;
    end
end


dbg = 1;
end


function sig_power = power_estimate(input_data)
sig_power = sum(input_data.*conj(input_data))/numel(input_data);
end

function [max_corr_val, max_corr_ind] = correlate(rx_signal, mask)

window_size = length(mask);
window_slide_dist = 1;
num_windows = ceil(((length(rx_signal) - window_size)/window_slide_dist) + 1);

corr_val = zeros(1, num_windows);

% autocorr in Taylor python code
for win = 1: num_windows
    
    start = (win-1)*window_slide_dist +1;
    fin = start + window_size - 1;
    
    corr_val(win) = abs(sum(rx_signal(start: fin)*mask'));
    
end
% figure()
% plot(corr_val)
% xlabel('Window index')
% ylabel('Correlation value')
% title(['Antenna ', num2str(ant), ' Correlation'])

[max_corr_val, max_corr_ind] = max(corr_val);

% if isempty(max_corr_ind)
%     error('not enough correlations, tx again or tx more data')
% end
%
% max_corr_ind = max(max_corr_ind);



dbg = 1;
end






