function [Fu,J] = fj_KdV_IM(ut,un,ABA,C,R,Mmm,Nm,ii,jj,k,p,order)
% Implicit midpoint time-integration
% %
Tut = tofu2(p,.5*(ut+un),order,R,Mmm,Nm,ii,jj);
dH = .5*C*(ut+un) - 1.5*Tut*(un+ut);
Fu = ut-un - k*ABA*dH;
%
I = speye(Mmm,Mmm);
if nargout > 1
   J = I - k*ABA*(.5*C - 3*Tut);
end