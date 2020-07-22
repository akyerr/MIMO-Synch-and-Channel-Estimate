function run_system(obj, snr_diag)


switch obj.sys_scenario
    case 'WIFIMIMOSM-A'
        obj.CP_type = 'Extended'; % 'Normal' or 'Extended'
        obj.num_ant = 2; % Tx and Rx
        obj.param_est = 'Ideal'; % 'Ideal' or 'Estimated
        obj.stream_size = obj.num_ant;
        
        if obj.diagnostic == 0
            obj.ebno_dB = [6, 7, 8, 9, 10, 14, 16, 20, 24];
        else
            obj.ebno_dB = snr_diag;
        end
        obj.num_subframes = [1, 10, 50, 100];
    case '4G5G'
        obj.CP_type = 'Normal'; % 'Normal' or 'Extended'
        obj.num_ant = 1; % Tx and Rx
        obj.param_est = 'Ideal'; % 'Ideal' or 'Estimated
        obj.stream_size = obj.num_ant;
        
        if obj.diagnostic == 0
            obj.ebno_dB = 30;
        else
            obj.ebno_dB = snr_diag;
        end
        obj.num_subframes = 1;
        
        
end
end

