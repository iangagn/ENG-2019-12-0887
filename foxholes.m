function [y] = foxholes(x)

a = [-32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32 -32 -16 0 16 32;
     -32 -32 -32 -32 -32 -16 -16 -16 -16 -16 0 0 0 0 0 16 16 16 16 16 32 32 32 32 32];    
    
nb_pts = size(x, 1);

y = zeros(nb_pts, 1);

for i = 1: nb_pts
    temp = 0;
    for k = 1:25
        temp = temp + 1/(k+(x(i,1) - a(1,k)).^6 + (x(i,2) - a(2,k)).^6);
    end
    y(i) = ((1/500) + temp).^(-1);
end

end