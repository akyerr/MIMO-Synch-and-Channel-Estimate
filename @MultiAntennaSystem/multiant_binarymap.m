function multiant_binarymap(obj, symb_pattern, binary_info)

% if strcmpi(obj.MIMO_method, 'MIMOTest')
    
    for ant = 1: obj.num_ant
        qpsk_data = lteSymbolModulate(binary_info(ant, :), obj.mod_type);
        
        symb_count = 0;
        datasymb_count = 0;
        
        for symb = 1: size(symb_pattern, 2)
            
            symb_count = symb_count + 1;
            
            if symb_pattern(ant, symb) == 0 % Zeros
                
                cmplx_data = zeros(1, length(obj.data_bin_ind));
                obj.tx_symbs(ant, obj.NFFT*(symb_count-1)+obj.data_bin_ind) = cmplx_data;
                
            elseif symb_pattern(ant, symb) == 1 % Synch
                
                cmplx_data = obj.ZChu(ant, :);
                obj.tx_symbs(ant, obj.NFFT*(symb_count-1)+obj.synch_bin_ind) = cmplx_data;
                
            elseif symb_pattern(ant, symb) == 2 % Data
                datasymb_count = datasymb_count + 1;
                qpsk_start = (datasymb_count-1)*obj.num_databins +1;
                qpsk_fin = qpsk_start + obj.num_databins - 1;
                cmplx_data = qpsk_data(qpsk_start: qpsk_fin);
                obj.tx_symbs(ant, obj.NFFT*(symb_count-1)+obj.data_bin_ind) = cmplx_data;
                
            end
            
        end
%                 figure();
%                 plot(real(obj.tx_symbs(ant, :)))
%                 dbg = 1;
    end
% else
%     error('Currently not supported')
% end

end
