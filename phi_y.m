function out = phi_y(x,order,node)
num = 0;
denom = 1;
for j = 1:order+1
    if j ~= node 
        numfac = 1;
        for k = 1:order+1
            if k ~= node && k~= j
                numfac = numfac*(x + 1 - 2*(k-1)/order);
            end
        end
        num = num + numfac;
        denom = denom*(2*(node-1)/order - 2*(j-1)/order);
    end
end
out = num/denom;
end