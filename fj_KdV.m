function [Fu,J] = fj_KdV(ut,un,ABA,C,R,Mmm,Nm,ii,jj,k,p,order)
% %
Tut = tofu2(p,ut,order,R,Mmm,Nm,ii,jj);
Tun = tofu2(p,un,order,R,Mmm,Nm,ii,jj);
dH = .5*C*(ut+un) - (Tun*un + .5*Tut*un + .5*Tun*ut + Tut*ut);
Fu = ut-un - k*ABA*dH;
%
I = speye(Mmm,Mmm);
if nargout > 1
   J = I - k*ABA*(.5*C - (Tun + 2*Tut));
end