function un = findun_mym(C,u,p,x,xo,H,Mm,R,Mmm,Nm,ii,jj,alg,order)
% Solving a minimization problem by the method of Lagrange multipliers to
% find an integral preserving transfer of u
% Finding us first from simple interpolation

us = interp1(xo,[u;u(1)],x(1:Mm),'pchip'); % methods: linear / spline / pchip

T = @(ut) tofu2(p,ut,order,R,Mmm,Nm,ii,jj);
Tus = T(us);

if alg == 1
    F = @(unl) [unl(1:Mm) - us - unl(Mm+1)*(C-3*Tus)*us; 
        hamkdv(p,unl(1:Mm),C,order) - H];
    options = optimset('Algorithm','trust-region-dogleg','Display','iter','MaxFunEvals',10000,'TolFun',1e-12);
elseif alg == 2
    F = @(unl) [unl(1:Mm) - us - unl(Mm+1)*(C-3*Tus)*us;
        hamkdv(p,unl(1:Mm),C,order) - H];
    options = optimset('Algorithm','trust-region-reflective','Display','iter','MaxFunEvals',10000,'TolFun',1e-12);
elseif alg == 9 % Feeding the Jacobian
    F = @(unl) fj_findun_mym(unl,us,Tus,Mm,p,C,H,R,Mmm,Nm,ii,jj,order);
    options = optimoptions(@fsolve,'Display','off',...
	'Algorithm','trust-region-dogleg',...
	'SpecifyObjectiveGradient',true,'PrecondBandWidth',0,'TolFun',1e-16);
else
    F = @(unl) [unl(1:Mm) - us - unl(Mm+1)*(C-3*Tus)*us; % 
        hamkdv(p,unl(1:Mm),C,order) - H];
    options = optimset('Display','off','TolFun',1e-10);
end
un = fsolve(F,[us;0],options);
un = un(1:Mm);