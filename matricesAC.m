function [A,C] = matricesAC(x,order)
E = doubint(order); % Integrate phi*phi
D = derivs(order); %Integrate phi'*phi'
M = order*(length(x)-1);
N = length(x)-1;

A_vec = zeros(1,N*length(E(:)));
C_vec = zeros(1,N*length(E(:)));
ntrips = 0;

l2g = repmat(1:order+1,N,1);
l2g = bsxfun(@plus,(0:order:(N-1)*order)',l2g);
l2g(end) = 1;
ii = repmat(l2g,1,order+1)'; % Replicate the array lg2 over N columns (then transpose)
ii = ii(:)'; % Make a single array of these indices, columnwise

l2g2 = repmat(1:order+1,order+1,1);
l2g2 = l2g2(:)';
l2g2 = repmat(l2g2,N,1);
l2g2 = bsxfun(@plus,(0:order:(N-1)*order)',l2g2);
l2g2(end,end-order:end) = 1;
l2g2 = l2g2';
jj = l2g2(:)';

for i = 1:N
   l2g = (order*(i-1)+1):(order*i+1);
   if i == N
      l2g(end) = 1;
   end
   L = length(l2g);
   
   kk = ntrips+1;
   kk2 = ntrips + L^2;
   A_vec(kk:kk2) = (x(i+1)-x(i))/2*E(:);
   C_vec(kk:kk2) = 2/(x(i+1)-x(i))*D(:);
   ntrips = ntrips + L^2;
end
A = sparse(ii,jj,A_vec,M,M);
C = sparse(ii,jj,C_vec,M,M);

end