function out = phi(x,order,node)
out = 1;
for i = 1:order+1
    if i ~= node 
        out = out*(.5*(x+1)*order - (i-1))/(node-i);
    end
end
end