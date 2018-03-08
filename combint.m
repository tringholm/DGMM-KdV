function G = combint(order)
% Integrate phi*phi'
% G = zeros(3,3);
G = zeros(order+1,order+1);
% phi_1 = @(y) y*(y-1)/2;
% phi_2 = @(y) -(y+1)*(y-1);
% phi_3 = @(y) y*(y+1)/2;
% phi_1y = @(y) y - 1/2;
% phi_2y = @(y) -2*y;
% phi_3y = @(y) y + 1/2;
% phi = {phi_1 phi_2 phi_3};
% phi_y = {phi_1y phi_2y phi_3y};
for i = 1:order+1
    for j = 1:order+1
%         f = @(y) phi{i}(y)*phi_y{j}(y);
        f = @(y) phi(y,order,i)*phi_y(y,order,j); % Order 2*order-1
        G(i,j) = gauss(f,[-1 1],order); % Exact for order 2*n-1 and less, so take n = order
    end
end