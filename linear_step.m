function y = linear_step(x)
 
% Description : LINEAR_STEP is a vectorized version of the step benchmark function.
%
% Usage       : LINEAR_STEP takes an m by n matrix as an input where m is the number of
%               points to evaluate and n is the number of dimensions of the space.
%
% Author      : Iannick Gagnon --> iannick.gagnon.1@ens.etsmtl.ca

y = 12 + sum(floor(x), 2);
      
end




    

 

