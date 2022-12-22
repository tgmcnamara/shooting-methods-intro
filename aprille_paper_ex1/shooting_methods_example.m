clear
clc

x0 = [0;0;0;0];
n = length(x0);
T = 1/60;
ode_opt = odeset('Jacobian',@circuit_jacobian);


run_manual = true;
if run_manual
% manual approach to see how many periods to run before it stabilizes
err_max = 10;
iters = 0;
tol = 1e-3;
max_errors = [];
max_iters = 1000;
% figure()
% hold
%while (err_max > tol) & (iters < max_iters)
for meep = 1:75
    iters = iters + 1;
    [t_sol, x_sol] = ode23t(@circuit_dynamics,[0,T], x0, ode_opt);
    errors = abs(x0-x_sol(end,:)')./max(abs(x0), [1;1;1;1]);
    err_max = norm(errors, inf);
    max_errors = [max_errors err_max];
    x0 = x_sol(end,:)';
    fprintf("Iter %d: error %.3f \n", iters, err_max)
%     plot(t_sol, x_sol(:,4))
end
disp(x0)
end

converged = false;
iters = 0;
max_iters = 100;
eps = 1e-4;
delta = 1e-4;
x0 = zeros(4,1);
x0_convergence = [-9.0734; 9.0555; .0090285; 9.1015];
k = 16667;
t_sol = linspace(0, T, k+1);
dt = t_sol(2)-t_sol(1);
while ~(converged)
    iters = iters + 1;
    norm_vec = max(abs(x0), 1);
%     x_sol = zeros(n, length(t));
%     x_sol(:,1)= x0;
%     % force constant step size
%     for i = 2:length(t)
%         [t_d,x_d] = ode23t(@circuit_dynamics,[t(i-1) t(i)], x_sol(:,i-1));
%         x_sol(:,i) = x_d(end,:);
%     end
    [t_sol2, x_sol2] = ode15s(@circuit_dynamics,[0,T], x0, ode_opt);
    x_sol = interp1(t_sol2, x_sol2, t_sol);
    x_sol = x_sol';
%     k = len(t_sol);
    xT = x_sol(:,end);
    cycle_error = norm((xT-x0)/norm_vec, inf);
    phi_inv = eye(n);
    phi = eye(n);
    for i = 2:length(t_sol)
%         index = k - i + 1;
%         fprintf("%d\n", i)
%         phi = (eye(n) - dt*circuit_jacobian(t_sol(index), x_sol(index)))\phi;
%         dt = t_sol(i)-t_sol(i-1);
        phi_inv = phi_inv * (eye(n) - dt*circuit_jacobian(t_sol(i), x_sol(:,i)));
    end
    % y = (eye(n)-phi)\(xT - phi*x0);
    y = (phi_inv - eye(n))\(phi_inv\xT-x0);
    iter_error = norm((y-x0)/norm_vec, inf);
    x0 = y;
    if cycle_error < eps && iter_error < delta
        converged = true;
        fprintf("Converged in %d iterations.\n", iters)
    elseif iters > max_iters
        fprintf("Max iterations reached. Failed to converge.\n")
        break
    else
        fprintf("Iter %d: error %.3e, %.3e\n", iters, cycle_error, iter_error)
    end
end