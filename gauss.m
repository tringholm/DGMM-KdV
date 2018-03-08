function out = gauss(f,x_bds,n)
[pts, wts] = lgwt(n,-1,1);
pts = pts(end:-1:1);
wts = wts(end:-1:1);
xpts = zeros(size(pts));
for i = 1:length(xpts)
   xpts(i) = (x_bds(2)-x_bds(1))/2*pts(i) + (x_bds(2)+x_bds(1))/2;
end
out = 0;
for i = 1:length(xpts)
   out = out + wts(i)*f(xpts(i)); 
end
out = out*(x_bds(2)-x_bds(1))/2;
end