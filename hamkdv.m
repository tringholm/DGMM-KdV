function H = hamkdv(x,u,C,order)
% Calculating the hamiltonian H = \int (.5*u'^2 - u^3)

R = tripint(order);
M = order*(length(x)-1);
N = length(x)-1;

T_vec = zeros(1,N*(order+1)^2);
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
   
   Ru = zeros(order+1,order+1); % u_k \sum phi_i phi_j phi_k
   for k = 1:order+1
        Ru = Ru + u(l2g(k))*R(:,:,k);
   end
   
   kk = ntrips+1;
   kk2 = ntrips + L^2;
   T_vec(kk:kk2) = (x(i+1)-x(i))/2*Ru(:); 
    ntrips = ntrips + L^2;
end
T = sparse(ii,jj,T_vec,M,M);
H = ((.5*C-T)*u)'*u;