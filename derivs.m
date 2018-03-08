function D = derivs(order)
% Integrate phi'*phi'
D = zeros(order+1,order+1);
% phi_1y = @(y) y - 1/2;
% phi_2y = @(y) -2*y;
% phi_3y = @(y) y + 1/2;
% phi = {phi_1y phi_2y phi_3y};
for i = 1:order+1
   for j = 1:order+1
       f = @(y) phi_y(y,order,i)*phi_y(y,order,j); % polys of order (order-1) -> total order of integrand is 2*order-2
       D(i,j) = gauss(f,[-1 1],order); % Exact for order 2*n-1 and less, so take n = order
   end
end
end