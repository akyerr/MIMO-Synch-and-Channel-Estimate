function multiant_binarymap(obj, symb_pattern, binary_info)


for ant = 1: obj.num_ant
    
    zerosymb_count = 0;
    synchsymb_count = 0;
    datasymb_count = 0;
    
    for symb = 1: size(symb_pattern, 2)
        
        if symb_pattern(ant, symb) == 0
            % Insert zeros
            % increase zero symb count
            zerosymb_count = zerosymb_count + 1;
        elseif symb_pattern(ant, symb) == 1
            % insert synch
            % increase synch symb count
            synchsymb_count = synchsymb_count + 1;
        elseif symb_pattern(ant, symb) == 2
            % QPSK map, insert data
            % increase data symbol count
            datasymb_count = datasymb_count + 1;
        end
    end
end


end
