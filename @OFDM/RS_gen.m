function RS = RS_gen(~, nRS)

if nargin == 2 % this means we want ref signals for all antenna
    rng('default')
    phases = [1, 3, 5, 7];
    RSall = zeros(nRS, 1);
    
    for i = 1: nRS
        idx = randperm(length(phases), 1);
        RSall(i) = exp(1i*2*pi/8*phases(idx));
    end
    RS = RSall;
end
