function [R,M,N,ii,jj] = tofuprep(x,order)
R = tripint(order);
M = order*(length(x)-1);
N = length(x)-1;

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