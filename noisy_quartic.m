function y = noisy_quartic(x)
 
% Description : NOISY_QUARTIC is a vectorized version of the quartic benchmark function optimized 
%               for low and high dimensional inputs.
%
% Usage       : NOISY_QUARTIC takes an m by n matrix as an input where m is the number of
%               points to evaluate and n is the number of dimensions of the space.
%
% Author      : Iannick Gagnon --> iannick.gagnon.1@ens.etsmtl.ca

    % Input dimensions
    [nb_input, nb_dim] = size(x);

    if nb_dim <= 50 % Low-dimensional implementation
        y = 0;
        for i = 1 : nb_dim
            y = y + i * x(:, i).^ 4;
        end
        y = y + rand();
    else
        y = sum(repmat(1:nb_dim, [nb_input, 1]) .* x.^4, 2) + rand(); % High-dimensional implementation
    end
    
end