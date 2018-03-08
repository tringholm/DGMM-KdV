function C = tripint(order)
% Integrate phi*phi*phi
C = zeros(order+1,order+1,order+1);
for i = 1:order+1
    for j = 1:order+1
        for k = 1:order+1
            f = @(y) phi(y,order,i)*phi(y,order,j)*phi(y,order,k);
            C(i,j,k) = gauss(f,[-1 1],2*order); % Exact for order 2*n-1 and less, so take n = 2*order
        end
    end
end
end