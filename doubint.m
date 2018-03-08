function E = doubint(order)
% Integrate phi*phi
% E = zeros(3,3);
E = zeros(order+1,order+1);
% phi_1 = @(y) y*(y-1)/2;
% phi_2 = @(y) -(y+1)*(y-1);
% phi_3 = @(y) y*(y+1)/2;
% phi = {phi_1 phi_2 phi_3};
for i = 1:order+1
    for j = 1:order+1
        %         f = @(y) phi{i}(y)*phi{j}(y);
        f = @(y) phi(y,order,i)*phi(y,order,j); % order 2*order
        E(i,j) = gauss(f,[-1 1],order+1); %Exact for order 2*n-1 and less, so take n = order + 1
    end
end
end