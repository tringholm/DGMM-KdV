function T = tofu2(x,u,order,R,M,N,ii,jj)
T_vec = zeros(1,N*(order+1)^2);
ntrips = 0;

for i = 1:N
   l2g = (order*(i-1)+1):(order*i+1);
   if i == N
      l2g(end) = 1;
   end
   L = length(l2g);
   
   Ru = zeros(order+1,order+1); 
   for k = 1:order+1
        Ru = Ru + u(l2g(k))*R(:,:,k);
   end
   
   kk = ntrips+1;
   kk2 = ntrips + L^2;
   T_vec(kk:kk2) = (x(i+1)-x(i))/2*Ru(:); 
   ntrips = ntrips + L^2;
end
T = sparse(ii,jj,T_vec,M,M);