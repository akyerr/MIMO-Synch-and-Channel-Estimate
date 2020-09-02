function vector = grid2vec(~, grid)
K = size(grid, 1); % number of frequency sub-carriers
L = size(grid, 2); % number of OFDM symbols
P = size(grid, 3); % Number of antennas

% vector will have dimensions L*K by P. L lots of K symbols appended to
% each other on each antenna.

vector = zeros(P, L*K);
for p = 1: P
    
    for l = 1: L
        start = (l-1)*K + 1;
        fin = start + (K-1);
        vector(p, start:fin) = grid(:, l, p);
    end
end

end

