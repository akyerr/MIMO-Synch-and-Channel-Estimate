function time_offset = synch_and_estimate(obj, synch_data, SNRdB)
SNR = 10^(SNRdB/10);
obj.stride_val = round(obj.CP/2);
obj.time_synch_ref = zeros(obj.num_ant,250,3);

for ant = 1: 1 % rx ant?
    obj.corr_obs = 0;
    chan_q = (reshape(obj.genie_channel_time(ant,1,:), [numel(obj.genie_channel_time(ant,1,:)),1]));
    obj.samp_delay = 4;
    obj.start_samp = obj.CP - obj.samp_delay;
    num_iter = ceil(size(obj.rx_waveform,2)/obj.stride_val);
    dlong = zeros(num_iter, 1);
    
    ptr_adj = 0;
    iter = 0; % P
    sym_count = 0;
    tap_delay = 5;
    x = zeros(1, tap_delay);
    rx_buffer_time = zeros(1, obj.NFFT*synch_data(1));
    while iter < num_iter
        iter = iter + 1;
        
        if obj.corr_obs == 0
            ptr_frame = (iter-1)*obj.stride_val + obj.start_samp + ptr_adj;
        elseif obj.corr_obs < 5
            ptr_frame =  ptr_frame + sum(synch_data)*(obj.CP+obj.NFFT);
        else
            ptr_frame  = round(XP(end,:)*b-obj.CP/4);
        end
        
        if  (synch_data(1)-1)*obj.samp_per_symb+obj.NFFT + ptr_frame -1 < ...
                size(obj.rx_waveform,2)
            for LL = 1: synch_data(1)
                start = (LL-1)*obj.samp_per_symb + ptr_frame;
                fin = (LL-1)*obj.samp_per_symb + ptr_frame+obj.NFFT-1;
                rx_buffer_time((LL-1)*obj.NFFT+1: LL*obj.NFFT) = obj.rx_waveform(ant, start: fin);
            end
            
            temp_vec1 = zeros(synch_data(1), obj.NFFT);
            
            for LL = 1: synch_data(1)
                temp_vec1(LL, 1: obj.NFFT) = fft(rx_buffer_time((LL-1)*obj.NFFT+1: LL*obj.NFFT), obj.NFFT);
            end
            
            synch00 = temp_vec1(:,obj.synch_bin_ind); % synchdat00
            synch0 = reshape(synch00.',[1,numel(synch00)]); % synchdat0
            Pest = sum(synch0.*conj(synch0))/length(synch0);
            synch = synch0/sqrt(Pest); % synchdat
            
            
            % delay estimation
            pmat0 = exp(1i*2*(pi/obj.NFFT)* (obj.synch_bin_ind-1).'*(0:obj.CP));
            pmat = repmat(pmat0,[synch_data(1), 1]);
            obj.del_mat=conj(obj.synch_ref)*(diag(synch)*pmat);
            [dmax, dmaxind0] = max( abs(obj.del_mat));
            
            dmaxind = dmaxind0 - 1;
            dlong(iter) = dmax;
            
            if dmax > 0.5*length(synch) || obj.corr_obs > 0 % if max correlation is more than halfway into the synch or if no corr obs found
                if dmaxind > round(0.75*obj.CP) % if max corr index is more than 75% of the way in the CP
                    
                    if obj.corr_obs == 0
                        ptr_adj = ptr_adj + round(0.5*obj.CP);
                        ptr_frame = (iter-1)*obj.stride_val + obj.start_samp + ptr_adj;
                    elseif obj.corr_obs < 5
                        ptr_frame = ptr_frame + round(0.5*obj.CP);
                    end
                    
                    for LL = 1: synch_data(1)
                        start = (LL-1)*obj.samp_per_symb + ptr_frame;
                        fin = (LL-1)*obj.samp_per_symb + ptr_frame+obj.NFFT-1;
                        rx_buffer_time((LL-1)*obj.NFFT+1: LL*obj.NFFT) = obj.rx_waveform(ant, start: fin);
                    end
                    
                    temp_vec1 = zeros(synch_data(1), obj.NFFT);
                    
                    for LL = 1: synch_data(1)
                        temp_vec1(LL, 1: obj.NFFT) = fft(rx_buffer_time((LL-1)*obj.NFFT+1: LL*obj.NFFT), obj.NFFT);
                    end
                    
                    synch00 = temp_vec1(:,obj.synch_bin_ind); % synchdat00
                    synch0 = reshape(synch00.',[1,numel(synch00)]); % synchdat0
                    Pest = sum(synch0.*conj(synch0))/length(synch0);
                    synch = synch0/sqrt(Pest); % synchdat
                    
                  
                    
                    % delay estimation
                    pmat0 = exp(1i*2*(pi/obj.NFFT)* (obj.synch_bin_ind-1).'*(0:obj.CP));
                    pmat = repmat(pmat0,[synch_data(1), 1]);
                    obj.del_mat=conj(obj.synch_ref)*(diag(synch)*pmat);
                    [dmax, dmaxind0] = max( abs(obj.del_mat));
                    
                    dmaxind = dmaxind0 - 1;
                    dlong(iter) = dmax;   
                end
                
                time_synch_ind = obj.time_synch_ref(ant, max(obj.corr_obs, 1), 1);
                
                if (ptr_frame - time_synch_ind) > (2*obj.CP + obj.NFFT+1) || obj.corr_obs == 0
                    obj.corr_obs = obj.corr_obs + 1;
                    obj.time_synch_ref(ant, obj.corr_obs, 1)= ptr_frame;
                    obj.time_synch_ref(ant, obj.corr_obs, 2)= dmaxind;
                    obj.time_synch_ref(ant, obj.corr_obs, 3)= dmax;
                    
                    ptr_synch0(rem(sym_count, tap_delay)+1) = sum(obj.time_synch_ref(ant, obj.corr_obs, 1:2)); %#ok
                    x(rem(sym_count, tap_delay)+1) = sym_count*sum(synch_data);
                    sym_count = sym_count + 1;
                    
                    x2 = x(1:min(obj.corr_obs, tap_delay));
                    xplus = [x2, sym_count*sum(synch_data)];
                    XP = [ones(length(xplus),1), xplus.'];
                    if obj.corr_obs > 3
                        y = (ptr_synch0(1: min(tap_delay,obj.corr_obs))).';
                        X = [ones(length(x2),1), x2.'];
                        b = X\y;
                      
                        dbg77=1;
                        
                    end
                    
                    
                    %%%%%Recovered Data with Delay Removed
                    data_recov = diag(synch)*pmat(:, dmaxind+1);
                    
                    hest1 = zeros(obj.NFFT, 1);
                    
                    tmp1 = data_recov.*(obj.synch_ref')./(1 + 1/SNR);
                    
                    LL = numel(tmp1)/synch_data(1);
                    
                    hest00 = reshape(tmp1, [LL, synch_data(1)]);
                    hest0 = hest00.';
                    hest = sum(hest0, 1)/synch_data(1); % averaging across adjacent synchs
                    
                    hest1(obj.synch_bin_ind) = hest;
                    
                    obj.est_chan_freqP(ant, obj.corr_obs, 1: length(hest1)) = hest1;
                    obj.est_chan_freqN(ant, obj.corr_obs, 1: length(hest)) = hest;
                    
                    
                    if obj.diagnostic == 1 && iter == 1
                        figure()
                        xax = (0:obj.NFFT - 1)*obj.fs/obj.NFFT;
                        yax1 = 20*log10(abs(hest1));
                        yax2 = 20*log10(abs(fft(chan_q, obj.NFFT)));
                        plot(xax, yax1, 'r', xax, yax2, 'b')
                    end
                    
                    
                end
                
            end
            
            
            
            
            
        end   
    end
end

time_offset = dmaxind;
dbg = 1;
end