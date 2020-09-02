function [RS_ind, indices0] = RS_indices(obj, num_symb_SF, index_type)
num_databins = obj.num_databins;
num_ant = obj.num_ant;


L = num_symb_SF/2;
l_init = [0, L-num_ant, L, 2*L-num_ant];

ports = 0: num_ant - 1;
indices0 = [];

for portno = ports
    l = l_init + portno;
    
    % spacing between 2 ref signals is 6 as per LTE
    temp_k = (0:6:num_databins - 3).' + [0, 3];
    len = length(temp_k);
    k = [temp_k, temp_k];
    k = k(:);
    l = repmat(l,len,1);
    l = l(:);
    p = repmat(portno,length(k),1);
    
    ind = [k, l, p];
    indices0 = [indices0; ind]; %#ok
    dbg = 1;
end

if (strcmpi(index_type, 'sub0'))
RS_ind = indices0;
elseif (strcmpi(index_type, 'sub1'))
    RS_ind = indices0 + 1;
elseif (strcmpi(index_type, 'lin'))
idx1 = indices0(:, 1)+1;
idx2 = indices0(:, 2)+1;
idx3 = indices0(:, 3)+1;

sz = [num_databins, num_symb_SF, num_ant];
RS_ind = sub2ind(sz, idx1, idx2, idx3);
end
end

