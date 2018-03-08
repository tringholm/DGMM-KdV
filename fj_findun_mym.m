function [F,J] = fj_findun_mym(unl,us,Tus,Mm,p,C,H,R,Mmm,Nm,ii,jj,order)

T = @(ut) tofu2(p,ut,order,R,Mmm,Nm,ii,jj);

F = [unl(1:Mm) - us - unl(Mm+1)*(C-3*Tus)*us;
    hamkdv(p,unl(1:Mm),C,order) - H]; 
%
if nargout > 1
   J = [speye(Mm),-(C-3*Tus)*us;
       ((C-3*T(unl(1:Mm)))*unl(1:Mm))',0];
end