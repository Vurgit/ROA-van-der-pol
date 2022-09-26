% Meng F ,  Wang D ,  Yang P , et al. Application of sum-of-squares method in
% estimation of region of attraction for nonlinear polynomial systems
% (August 2019)[J]. IEEE Access, 2020, PP(99):1-1.
clear;
pvar x1 x2
pvar beta gamma
x = [x1; x2];
% prog = sosprogram(x);
f = [-x2; x1 + (x1^2 - 1) * x2];
jac = jacobian(f, x);
degV = 4;
V0 = 0.65217 * x1^2 - 0.43478 * x1 * x2 + 0.43478 * x2^2;
p0 = x.' * x;
monos = monomials(x, 1 : degV);
% [prog, V] = sossosvar(prog, monos);
% dotV = jacobian(V, vars) * f;
% dotV = sospolyvar(prog, monos, 'wscoeff');
L1 = 1e-5 * (x).' * x;
L2 = L1;

% [prog, s1] = sossosvar(prog, monos);
% [prog, s1] = sossosvar(prog, monos);
% [prog, s2] = sossosvar(prog, monos);
V = V0;
p = p0;
gamma0 = 1;
beta0 = 1;
gamma_sol = gamma0;
beta_sol = beta0;
iter = 0;
while 1
    
    beta_sol_old = beta_sol;
    
    % solve for gamma and s2 when V fixed
    dotV = jacobian(V, x) * f;
    prog1_1 = sosprogram(x);
    [prog1_1, s2] = sossosvar(prog1_1, monos);
    prog1_1 = sosineq(prog1_1, -(dotV + L2 + s2 * (gamma_sol - V)) );
    prog1_1 = sossolve(prog1_1);
    s2_sol = sosgetsol(prog1_1, s2);
    
    prog1_2 = sosprogram(x, gamma);
    prog1_2 = sosineq(prog1_2, -(dotV + L2 + s2_sol * (gamma - V)) );
    prog1_2 = sossetobj(prog1_2, -gamma);
    prog1_2 = sossolve(prog1_2);
    gamma_sol = sosgetsol(prog1_2, gamma);
    
    
    
    prog2_1 = sosprogram(x);
    [prog2_1, s1] = sossosvar(prog2_1, monos);
    prog2_1 = sosineq(prog2_1, -((V - gamma_sol) + s1 * (beta_sol - p)) );
    prog2_1 = sossolve(prog2_1);
    s1_sol = sosgetsol(prog2_1, s1);
    
    prog2_2 = sosprogram(x, beta);
    prog2_2 = sosineq(prog2_2, -((V - gamma_sol) + s1_sol * (beta - p)) );
    prog2_2 = sossetobj(prog2_2, -beta);
    prog2_2 = sossolve(prog2_2);
    beta_sol = sosgetsol(prog2_2, beta);
    if abs( double(beta_sol - beta_sol_old) ) < 1e-3
        break;
    end
    iter = iter + 1;
    
    beta_sol
    
    
    prog3 = sosprogram(x);
    [prog3, V] = sossosvar(prog3, monos);
    dotV = jacobian(V, x) * f;
    %     prog3 = sosineq(V);
    prog3 = sosineq(prog3, -((V - gamma_sol) + s1_sol * (beta_sol - p)) );
    prog3 = sosineq(prog3, -(dotV + L2 + s2_sol * (gamma_sol - V)) );
    prog3 = sossolve(prog3);
    V = sosgetsol(prog3, V);
    V = V / gamma_sol;
end

domain = [-3 3 -3 3];
[C,ph(1)]=pcontour(V,double(gamma_sol),domain,'g');
hold on

