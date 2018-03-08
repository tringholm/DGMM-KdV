function [un,p,H] = KdV_AVF(M,k,Tmax,L,method,order,moving,monitor,doplot,u_analytic)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IF USING THIS CODE FOR RESEARCH PURPOSES, PLEASE CITE OUR ARTICLE     %
% Eidnes, S., Owren, B. & Ringholm, T. Adv Comput Math (2017).          %
% https://doi.org/10.1007/s10444-017-9562-8                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the 1D Korteweg?de Vries equation numerically using moving mesh and the AVF scheme for time stepping.
%
% Initial solution: u_analytic(x,0)
% Periodic boundary conditions: u(L,t) = u(-L,t)
% Preserving H = \int(.5*u'^2 - u^3)dx
% Using Legendre-Gauss quadrature to evaluate integrals, employing the code lgwt.m written by Greg von Winckel
%
% Input: 
%
% M             - Number of spatial discretization points (-1)
% k             - Time step
% Tmax          - End time
% L             - Interval [-L,L]
% method        - 1 for AVF DG method, 2 for implicit midpoint method
% order         - order of the basis functions
% moving        - true if using moving mesh
% monitor       - 1 for arc-length monitor function, 2 for curvature
% doplot        - true if producing plots while running
% u_analytic    - analytic solution for u, for finding the initial value and for testing purposes
%
% Output: 
%
% un            - solution u at time T
% p             - spatial discretization points at time T
% H             - Hamiltonian values at all time steps

p = linspace(-L,L,M+1)';
P = 2*L;
N = Tmax/k;

up = u_analytic(p,0); % Initial solution
if moving % Creating the initial mesh
    fa = 100; % k^2 in the generalized monitor function
    lv = zeros(1,M);
    e2 = ones(M,1);
    Hip = spdiags([-e2 e2], 0:1, M, M+1);
    hip = Hip*p;
    if monitor == 1
        lv(1) = sqrt(hip(1)^2+fa*(up(2)-up(1))^2);
        j = 2:M-1;
        lv(j) = sqrt(hip(j).^2+fa*(up(j+1)-up(j)).^2);
        lv(M) = sqrt(hip(M)^2+fa*(up(1)-up(M))^2);
    elseif monitor == 2
        lv(1) = .5*(hip(1)+hip(2)).*nthroot(1+fa*(2*(hip(1).*up(2)-(hip(1)+hip(2)).*up(1))./(hip(1).*hip(2).^2+hip(2).*hip(1).^2)).^2,4);
        j = 2:M-2;
        ddux = 2*(hip(j).*up(j+1)-(hip(j)+hip(j+1)).*up(j)+hip(j+1).*up(j-1))./(hip(j).*hip(j+1).^2+hip(j+1).*hip(j).^2); % 1st order approximation
        lv(j) = .5*(hip(j)+hip(j+1)).*nthroot(1+fa*ddux.^2,4);
        lv(M-1) = .5*(hip(M-1)+hip(M)).*nthroot(1+fa*2*(-(hip(M-1)+hip(M)).*up(M-1)+hip(M).*up(M-2))./(hip(M-1).*hip(M).^2+hip(M).*hip(M-1).^2).^2,4);
        lv(M) = .5*(hip(M)+hip(1)).*nthroot(1+fa*(2*(hip(M).*up(1)+hip(1).*up(M-1))./(hip(M).*hip(1).^2+hip(1).*hip(M).^2)).^2,4);
    end
    Le = sum(lv);
    l = Le/M;
    po = p;
    sn = 0;
    n = 1;
    for j = 1:M-1
        while sn + lv(n) < l*j
            sn = sn + lv(n);
            n = n+1;
        end
        if n == M
            p(j+1) = po(M) + (l*j-sn)*(po(1)+P-po(M))/lv(M);
        else
            p(j+1) = po(n) + (l*j-sn)*(po(n+1)-po(n))/lv(n);
        end              
    end
end

Mm = order*M;
x = zeros(Mm + 1,1); % Contains all points, including virtual nodes inside elements
for i = 1:M
    for j = 1:order
        x(order*(i-1) + j) = p(i) + (j-1)*(p(i+1)-p(i))/order;
    end
end
x(end) = p(end);

un = u_analytic(x(1:Mm),0);
mh = max(un);

[~,C] = matricesAC(p,order);
H0 = hamkdv(p,un,C,order);
H = [H0,zeros(1,N)];

if method == 1 % AVF DG method
    [R,Mmm,Nm,ii,jj] = tofuprep(p,order);
    if doplot
        plot(x,[un;un(1)],'b','LineWidth',1)
        axis([-100 100 -.5 3.5])
        hold on
        plot(x,[un;un(1)],'r','LineWidth',1)
        hold off
        pause(0.01)
    end
    for i = 1:N
        if moving
            hip = Hip*p;
            up = un(1:order:end);
            if monitor == 1
                lv(1) = sqrt(hip(1)^2+fa*(up(2)-up(1))^2);
                j = 2:M-1;
                lv(j) = sqrt(hip(j).^2+fa*(up(j+1)-up(j)).^2);
                lv(M) = sqrt(hip(M)^2+fa*(up(1)-up(M))^2);
            elseif monitor == 2
                lv(1) = .5*(hip(1)+hip(2)).*nthroot(1+fa*(2*(hip(1).*up(2)-(hip(1)+hip(2)).*up(1))./(hip(1).*hip(2).^2+hip(2).*hip(1).^2)).^2,4);
                j = 2:M-2;
                ddux = 2*(hip(j).*up(j+1)-(hip(j)+hip(j+1)).*up(j)+hip(j+1).*up(j-1))./(hip(j).*hip(j+1).^2+hip(j+1).*hip(j).^2); % 1st order approximation
                lv(j) = .5*(hip(j)+hip(j+1)).*nthroot(1+fa*ddux.^2,4);
                lv(M-1) = .5*(hip(M-1)+hip(M)).*nthroot(1+fa*2*(-(hip(M-1)+hip(M)).*up(M-1)+hip(M).*up(M-2))./(hip(M-1).*hip(M).^2+hip(M).*hip(M-1).^2).^2,4);
                lv(M) = .5*(hip(M)+hip(1)).*nthroot(1+fa*(2*(hip(M).*up(1)+hip(1).*up(M-1))./(hip(M).*hip(1).^2+hip(1).*hip(M).^2)).^2,4);
            end
            Le = sum(lv);
            l = Le/M;
            xo = x;
            po = p;
            sn = 0;
            n = 1;
            for j = 1:M-1
                while sn + lv(n) < l*j
                    sn = sn + lv(n);
                    n = n+1;
                end
                if n == M
                    p(j+1) = po(M) + (l*j-sn)*(po(1)+P-po(M))/lv(M);
                else
                    p(j+1) = po(n) + (l*j-sn)*(po(n+1)-po(n))/lv(n);
                end              
            end
        end
        for it = 1:M
            for j = 1:order
                x(order*(it-1) + j) = p(it) + (j-1)*(p(it+1)-p(it))/order;
            end
        end
        x(end) = p(end);
%        
        [A,C] = matricesAC(p,order);
        if moving
            alg = 9; 
            un = findun_mym(C,un,p,x,xo,H0,Mm,R,Mmm,Nm,ii,jj,alg,order); 
        end
        B = matrixD(p,order);
        ABA = (A^-1)*B*(A^-1);
        Fu = @(ut) fj_KdV(ut,un,ABA,C,R,Mmm,Nm,ii,jj,k,p,order);
        options = optimoptions(@fsolve,'Display','off',...
        	'Algorithm','trust-region-dogleg',...
            'SpecifyObjectiveGradient',true,'PrecondBandWidth',0,'TolFun',1e-12);
        un = fsolve(Fu,un,options);
        H(i+1) = hamkdv(p,un,C,order);
        if doplot
            plot(x,[un;un(1)],'b','LineWidth',1)
            hold on
            plot(x,u_analytic(x,i*k),'r','LineWidth',1)
            axis([-L L -.5 mh+.5])
            hold off
            pause(0.01)
        end
    end
elseif method == 2 % Implicit midpoint
    [R,Mmm,Nm,ii,jj] = tofuprep(p,order);
    if doplot
        plot(x,[un;un(1)],'b','LineWidth',1)
        axis([-100 100 -.5 3.5])
        hold on
        plot(x,[un;un(1)],'r','LineWidth',1)
        hold off
        pause(0.01)
    end
    for i = 1:N
        if moving
            hip = Hip*p;
            up = un(1:order:end);
            if monitor == 1
                lv(1) = sqrt(hip(1)^2+fa*(up(2)-up(1))^2);
                j = 2:M-1;
                lv(j) = sqrt(hip(j).^2+fa*(up(j+1)-up(j)).^2);
                lv(M) = sqrt(hip(M)^2+fa*(up(1)-up(M))^2);
            elseif monitor == 2
                lv(1) = .5*(hip(1)+hip(2)).*nthroot(1+fa*(2*(hip(1).*up(2)-(hip(1)+hip(2)).*up(1))./(hip(1).*hip(2).^2+hip(2).*hip(1).^2)).^2,4);
                j = 2:M-2;
                ddux = 2*(hip(j).*up(j+1)-(hip(j)+hip(j+1)).*up(j)+hip(j+1).*up(j-1))./(hip(j).*hip(j+1).^2+hip(j+1).*hip(j).^2); % 1st order approximation
                lv(j) = .5*(hip(j)+hip(j+1)).*nthroot(1+fa*ddux.^2,4);
                lv(M-1) = .5*(hip(M-1)+hip(M)).*nthroot(1+fa*2*(-(hip(M-1)+hip(M)).*up(M-1)+hip(M).*up(M-2))./(hip(M-1).*hip(M).^2+hip(M).*hip(M-1).^2).^2,4);
                lv(M) = .5*(hip(M)+hip(1)).*nthroot(1+fa*(2*(hip(M).*up(1)+hip(1).*up(M-1))./(hip(M).*hip(1).^2+hip(1).*hip(M).^2)).^2,4);
            end
            Le = sum(lv);
            l = Le/M;
            xo = x;
            po = p;
            sn = 0;
            n = 1;
            for j = 1:M-1
                while sn + lv(n) < l*j
                    sn = sn + lv(n);
                    n = n+1;
                end
                if n == M
                    p(j+1) = po(M) + (l*j-sn)*(po(1)+P-po(M))/lv(M);
                else
                    p(j+1) = po(n) + (l*j-sn)*(po(n+1)-po(n))/lv(n);
                end              
            end
        end
        for it = 1:M
            for j = 1:order
                x(order*(it-1) + j) = p(it) + (j-1)*(p(it+1)-p(it))/order;
            end
        end
        x(end) = p(end);
%   
        [A,C] = matricesAC(p,order);
        if moving % if moving mesh
            un = interp1(xo,[un;un(1)],x(1:Mm),'pchip'); % methods: linear / spline / pchip
        end
        B = matrixD(p,order);
        ABA = (A^-1)*B*(A^-1);
        Fu = @(ut) fj_KdV_IM(ut,un,ABA,C,R,Mmm,Nm,ii,jj,k,p,order);
        options = optimoptions(@fsolve,'Display','off',...
        	'Algorithm','trust-region-dogleg',...
            'SpecifyObjectiveGradient',true,'PrecondBandWidth',0,'TolFun',1e-12);
        un = fsolve(Fu,un,options);
        H(i+1) = hamkdv(p,un,C,order);
        if doplot
            plot(x,[un;un(1)],'b','LineWidth',1)
            hold on
            plot(x,u_analytic(x,i*k),'r','LineWidth',1)
            axis([-L L -.5 mh+.5])
            hold off
            pause(0.01)
        end
    end
else
    disp('Error: Wrong method input')
end
end