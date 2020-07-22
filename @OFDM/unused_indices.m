function unused_ind = unused_indices(obj, RS_ind)
num_ant = obj.num_ant;
num_databins = obj.num_databins;
% Separate reference signal indices on on each antenna. Each column has the
% linear indices for an antenna
rs_ant = reshape(RS_ind, [], num_ant);

unused_indices = cell(num_ant, 1);
offsets = num_databins: num_databins: num_ant*num_databins;
for ant = 1: num_ant
    offsets = offsets - num_databins;
    offsets_without_0 = offsets;
    offsets_without_0(offsets_without_0 == 0) = [];
    unused_bin_offsets = offsets_without_0;
    
    unused_indices{ant} = rs_ant(:, ant) + unused_bin_offsets;
end

unused_ind = reshape(cell2mat(unused_indices), numel(cell2mat(unused_indices)), 1);
dbg = 1;
end

%



