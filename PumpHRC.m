% PumpHRC.m for IFAC WC
% by Muhammad Hilmi, ITK NTNU
% at 23-11-2025

yalmip('clear');
clear all; close all;

rng(100);
load('FeedDist.mat');

% Model parameters
C = 1*eye(3); dt = 0.1; 
J = 5*eye(3); H = 0.1*eye(3);
K = 0.01*eye(3); Phi = 1*eye(3);
Ja = 5*eye(6); Ha = [zeros(3) 0.1*eye(3)];
Ca = [eye(3) zeros(3)]; Pa = [zeros(3) eye(3); zeros(3,6)];

% RTO and MPC configuration
Np = 2; Nq = 100; Ns = 5; % Change for faster/slower simulation!!
nx = 3; nu = 3; ny = 3; nc = 2; 

% State and input bounds
umin = [1600; 1500; 0]; umax = [3800; 4700; 1]; 
xmin = [0; 0; 0]; xmax = [276.2; 551.7; 1400];

% System parameters
pe = 206.2; lambda_eps = 1e3; qnom = 458; hm = 3; fp = 3;
pa = 32.02; pb = 217.45; rho = 1037; g = 9.81; rm = 0.0274;

% YALMIP variables
u_s = sdpvar(nu, Ns+1, 'full'); r = sdpvar(ny, 1, 'full');
x_s = sdpvar(nx, Ns+1, 'full'); xbar_s = sdpvar(nx,1,'full');
eps_q = sdpvar(1, Ns, 'full'); eps_p = sdpvar(1, Ns, 'full');
u = sdpvar(nu, Np,'full'); ubar = sdpvar(nu, 1,'full');
x = sdpvar(nx, Np+1, Nq, 'full'); xbar = sdpvar(nx, 1, 'full');
d_param = sdpvar(nx, Nq, 'full'); d = sdpvar(nx,1,'full');
X0_param = sdpvar(nx, Nq, 'full'); alpha = binvar(1, 1, 'full');
pd = sdpvar(1, 1, 'full'); qd = sdpvar(1, 1, 'full'); 

% RTO formulation
obj = 0;
cons = [];
M3 = umax(3);

% Receding horizon
for k = 1:Ns
    % Dynamic constraint
    cons = [cons, ...
        x_s(:,k+1) == cascade_approx(x_s(:,k), u_s(:,k), pd, qd), ...
        r == C * xbar_s];

    % Optimization cost
    obj = obj + (x_s(:,k) - xbar_s)' * (x_s(:,k) - xbar_s) + ...
        lambda_eps * (eps_p(k)^2 + eps_q(k)^2);
    if k < Ns
        obj = obj + (u_s(:,k+1) - u_s(:,k))' * (u_s(:,k+1) - u_s(:,k));
    end

    % Safety constraint
    h_npsh = pa + pb * ((qd - x_s(3,k))/qnom)^2; % NPSHr approximation
    p_npsh = rho * g * (h_npsh + hm) * 1e-5; dp = rm * qd + fp;
    cons = [cons, p_npsh - (x_s(1,k) - dp) <= eps_p(k)];
    cons = [cons, x_s(3,k) - qd/2 <= eps_q(k)];
    cons = [cons, eps_p(k) >= 0]; cons = [cons, eps_q(k) >= 0];

    % State and input constraints
    cons = [cons, xmin <= x_s(:,k) <= xmax];
    cons = [cons, umin <= u_s(:,k) <= umax];
    cons = [cons, xmin <= xbar_s <= xmax];

    % Mixed-integer constraint
    cons = [cons, u_s(3,k) <= M3 * alpha];
    cons = [cons, implies(alpha, x_s(1,k) >= pe)];
    cons = [cons, implies(~alpha, x_s(1,k) <= pe - 1e-3)];
end

% MPC formulation
objective = 0;
constraints = [];

% Offset-free constraints
constraints = [constraints, ...
    xbar == cascade_approx(xbar, ubar, pd, qd) + Phi * d * dt];
constraints = [constraints, umin <= ubar <= umax];

% Scenario and prediction horizon
for s = 1:Nq
    % Initial state
    constraints = [constraints, x(:,1,s) == X0_param(:,s)];
    for k = 1:Np
        % Expected cost
        objective = objective + 1/Nq * ...
            ((x(:,k,s) - xbar)' * (x(:,k,s) - xbar) + ...
            (u(:,k) - ubar)' * (u(:,k) - ubar));

        % Dynamic constraint
        constraints = [constraints, ...
            x(:,k+1,s) == cascade_approx(x(:,k,s), u(:,k), pd, qd) + ...
                Phi * d_param(:,s) * dt];

        % State and input constraints
        constraints = [constraints, xmin <= x(:,k,s) <= xmax];
        constraints = [constraints, umin <= u(:,k) <= umax];

        % Mixed-integer constraint
        constraints = [constraints, u(3,k) <= M3 * alpha];
        constraints = [constraints, implies(alpha, x(1,k,s) >= pe)];
        constraints = [constraints, implies(~alpha, x(1,k,s) <= pe - 1e-3)];
    end
end

% Finalize optimizer
param_in = {x_s(:,1), u_s(:,1), r, pd, qd};
sol_out = {xbar_s, u_s}; opt = sdpsettings('solver','gurobi');
rto = optimizer(cons, obj, opt, param_in, sol_out);

% Finalize controller
parameters_in = {X0_param, xbar, d_param, d, pd, qd};
solutions_out = {u, ubar}; options = sdpsettings('solver','gurobi');
controller = optimizer(constraints, objective, options, ...
    parameters_in, solutions_out);

% Initial conditions
t = 0:dt:300; xact = [113; 220; 0]; tact = [0; 0; 0];
yact = xact; xhat = xact; that = tact; xhat2 = [xact; tact];
xactprev = xact; xhatprev = xact; uact = [3063; 3514; 0];
env = [linspace(0,1,1200), ones(1,601), linspace(1,0,1200)];
xactArray = []; tactArray = []; yactArray = []; uactArray = [];
xhatArray = []; thatArray = []; xhat2Array = []; 

% Uncertainty initialization
xhat_samples = zeros(3, Nq); that_samples = zeros(3, Nq);
xhat_mean = zeros(3, length(t)); that_mean = zeros(3, length(t)); 
xhat_MC = zeros(3, Nq, length(t)); that_MC = zeros(3, Nq, length(t));
xsprev = xact; xhatprev_samples = xhat_samples;

tic;
hbar = waitbar(0, "Simulation progress...");

% Simulation
for i = 1:length(t)
    % Noise signals
    vact = -1 + 2 * sum(rand(3,100),2) / 100; wact = vact;

    % Fault signals
    tact = [0.1 * env(i) * sin(0.03 * pi * t(i)); ...
        0.2 * env(i) * cos(0.02 * pi * t(i)); ...
        0.25 * (1 - exp(-t(i) / 60))];

    % Supervisory RTO
    ins = {xhat, uact, ract(:,i), pf(i), qf(i)};
    [outs, diags] = rto{ins}; xp = outs{1};
    if diags ~= 0; xss = ract(:,i); else; xss = xp; end

    % Stochastic offset-free NMPC
    inputs = {xhat_samples, xss, that_samples, that, pf(i), qf(i)};
    [outputs, diagnostics] = controller{inputs}; up = outputs{1};
    if diagnostics ~= 0; uact = uact; else; uact = up(:,1); end

    % Store variables
    xactArray = [xactArray xact]; tactArray = [tactArray tact];
    xhatArray = [xhatArray xhat]; thatArray = [thatArray that];
    yactArray = [yactArray yact]; uactArray = [uactArray uact];
    xhat2Array = [xhat2Array xhat2];

    % Dynamic evolution
    fact = cascade(xact, uact, pf(i), qf(i));
    xact = fact + Phi * tact * dt + vact * dt;
    yact = C * xact + wact;
   
    % Nonlinear observer (joint)
    G = -inv(J) * H' * C;
    P = 2 * inv(Phi * Phi + K);
    fhat = cascade(xhat, uact, pf(i), qf(i));
    xhat = fhat + Phi * that * dt + inv(J) * H' * (yact - C * xhat) * dt;
    that = that + Phi' * P * ...
        ((xact - xhat) - K * (xactprev - xhatprev) - (fact - fhat));

    % Nonlinear observer (conventional)
    fhat2 = cascade(xhat2, uact, pf(i), qf(i));
    xhat2 = fhat2 + Pa * xhat2 * dt + inv(Ja) * Ha' * (yact - Ca * xhat2) * dt;

    % Monte-Carlo uncertainty propagation
    for s = 1:Nq
        vs = -1 + 2 * sum(rand(3,100),2) / 100; ws = vs;
        xs = fact + Phi * tact * dt + vs * dt;
        ys = C * xs + ws;

        % Update for sampled observer
        fhat_MC = cascade(xhat_samples(:,s), uact, pf(i), qf(i));

        % State estimate propagation
        xhat_samples(:,s) = fhat_MC + ...
              Phi * that_samples(:,s) * dt + ...
              inv(J) * H' * (ys - C * xhat_samples(:,s)) * dt;

        % Fault estimate propagation
        that_samples(:,s) = that_samples(:,s) + ...
              Phi' * P * ((xs - xhat_samples(:,s)) ...
                       - K * (xsprev - xhatprev_samples(:,s)) ...
                       - (fact - fhat_MC));

        % Store for uncertainty
        xhat_MC(:,s,i) = xhat_samples(:,s);
        that_MC(:,s,i) = that_samples(:,s);

        xsprev = xact;
    end

    % Mean and Standard deviation over samples
    xhat_mean(:,i) = mean(xhat_samples, 2);
    that_mean(:,i) = mean(that_samples, 2);
    xhat_std(:,i) = std(xhat_samples, 0, 2);
    that_std(:,i) = std(that_samples, 0, 2);

    xactprev = xact; xhatprev = xhat;
    xhatprev_samples = xhat_samples;

    % Norm errors
    nxhat(i) = norm(xact - xhat); nxhat2(i) = norm(xact - xhat2(1:3)); 
    nthat(i) = norm(tact - that); nthat2(i) = norm(tact - xhat2(4:end));

    waitbar(i * dt / t(end), hbar);
end

close(hbar);
T = toc;

%% Visualization
% State figure
fh = figure(1);
fh.Position = [220 300 620 470];

subplot(3,1,1);
plot(t, xactArray(1,:), 'LineWidth', 10, 'Color', "#466C95"); hold on;
plot(t, xhatArray(1,:), ':', 'LineWidth', 5, 'Color', "#5DAE8B");
legend('State Trajectory', 'State Estimation', 'Location', 'best');
ylabel('$p_f\;(\mathrm{bar})$', 'Interpreter','latex');
set(gca, 'FontSize', 18); hold off; 

subplot(3,1,2);
plot(t, xactArray(2,:), 'LineWidth', 10, 'Color', "#466C95"); hold on;
plot(t, xhatArray(2,:), ':', 'LineWidth', 5, 'Color', "#5DAE8B");
ylabel('$p_i\;(\mathrm{bar})$', 'Interpreter','latex');
set(gca, 'FontSize', 18); hold off;  

subplot(3,1,3);
plot(t, xactArray(3,:), 'LineWidth', 10, 'Color', "#466C95"); hold on;
plot(t, xhatArray(3,:), ':', 'LineWidth', 5, 'Color', "#5DAE8B");
xlabel('$t\;(\mathrm{s})$', 'Interpreter','latex');
ylabel('$q_r\;(\mathrm{m^3/h})$', 'Interpreter','latex');
set(gca, 'FontSize', 18); hold off; 

% Fault figure
fh = figure(2);
fh.Position = [860 520 620 470];

% Uncertainty envelope
tunc_ub = thatArray + 5 * that_std;
tunc_lb = thatArray - 5 * that_std;

subplot(3,1,1);
plot(t, tactArray(1,:), 'LineWidth', 10, 'Color', "#DA5C53"); hold on;
plot(t, thatArray(1,:), ':', 'LineWidth', 5, 'Color', "#5DAE8B");
fill([t fliplr(t)], [tunc_ub(1,:) fliplr(tunc_lb(1,:))], ...
    hex2rgb("#5DAE8B"), 'FaceAlpha', 0.3, 'EdgeColor','none');
ylabel('$\theta_{p_f}\;(\mathrm{bar})$', 'Interpreter','latex');
set(gca, 'FontSize', 18); ylim([-0.2 0.2]); hold off;

subplot(3,1,2);
plot(t, tactArray(2,:), 'LineWidth', 10, 'Color', "#DA5C53"); hold on;
plot(t, thatArray(2,:), ':', 'LineWidth', 5, 'Color', "#5DAE8B");
fill([t fliplr(t)], [tunc_ub(2,:) fliplr(tunc_lb(2,:))], ...
    hex2rgb("#5DAE8B"), 'FaceAlpha', 0.3, 'EdgeColor','none');
ylabel('$\theta_{p_i}\;(\mathrm{bar})$', 'Interpreter','latex');
set(gca, 'FontSize', 18); ylim([-0.25 0.25]); hold off;

subplot(3,1,3);
plot(t, tactArray(3,:), 'LineWidth', 10, 'Color', "#DA5C53"); hold on;
plot(t, thatArray(3,:), ':', 'LineWidth', 5, 'Color', "#5DAE8B");
fill([t fliplr(t)], [tunc_ub(3,:) fliplr(tunc_lb(3,:))], ...
    hex2rgb("#5DAE8B"), 'FaceAlpha', 0.3, 'EdgeColor','none');
legend('Fault Signal', 'Estimation Signal', ...
    'Orientation', 'horizontal', 'Location', 'best');
xlabel('$t\;(\mathrm{s})$', 'Interpreter','latex')
ylabel('$\theta_{q_r}\;(\mathrm{m^3/h})$', 'Interpreter','latex');
set(gca, 'FontSize', 18); ylim([-0.2 0.6]); hold off;

% Comparison figure
fh = figure(3);
fh.Position = [850 80 560 320];

subplot(2,1,1);
plot(t, nxhat, ':', 'LineWidth', 10, 'Color', "#5DAE8B"); hold on;
plot(t, nxhat2, ':', 'LineWidth', 10, 'Color', "#FEB21A");
ylabel('$\|\hat{x}\|$', 'Interpreter','latex');
set(gca, 'FontSize', 18); hold off;

subplot(2,1,2);
plot(t, nthat, ':', 'LineWidth', 8, 'Color', "#5DAE8B"); hold on;
plot(t, nthat2, ':', 'LineWidth', 8, 'Color', "#FEB21A");
legend('Proposed Estimation', 'Conventional Estimation', ...
    'Orientation', 'horizontal', 'Location', 'best');
xlabel('$t\;(\mathrm{s})$', 'Interpreter','latex')
ylabel('$\|\hat{\theta}\|$', 'Interpreter','latex');
set(gca, 'FontSize', 18); hold off;

% Cacade pump system model
function xk = cascade(xk, uk, pd, qd)
    rho = 1037; g = 9.81; rm = 0.0274; fp = 3;
    pe = 206.2; cv = 11.95;

    % Pump head coefficients
    pf1 = 8.45e-2; pf2 = 2.48e-4; pf3 = -7.91e-4; 
    pi1 = 0.3075; pi2 = 1.76e-3; pi3 = -4.63e-3;

    hf = pump(uk(1), qd, pf1, pf2, pf3);
    hi = pump(uk(2), qd-xk(3), pi1, pi2, pi3);

    pci = max(xk(1) - pe, 0);
    dp = rm * qd + fp;

    xk(1) = pd + rho * g * hf * 1e-5;
    xk(2) = (xk(1) - dp) + rho * g * hi * 1e-5;
    xk(3) = cv * uk(3) * (pci * 1e5 / rho)^0.5;
end

% Local linearized cascade pump system
function xk = cascade_approx(xk, uk, pd, qd)
    rho = 1037; g = 9.81; rm = 0.0274; fp = 3;
    pe = 206.2; cv = 11.95; pnom = 244.6;

    % Pump head coefficients
    pf1 = 8.45e-2; pf2 = 2.48e-4; pf3 = -7.91e-4; 
    pi1 = 0.3075; pi2 = 1.76e-3; pi3 = -4.63e-3;

    hf = pump(uk(1), qd, pf1, pf2, pf3);
    hi = pump(uk(2), qd-xk(3), pi1, pi2, pi3);

    alpha = 1e5 / (2 * rho * sqrt((pnom - pe) * 1e5 / rho));
    pci = max(xk(1) - pnom, 0);
    dp = rm * qd + fp;

    xk(1) = pd + rho * g * hf * 1e-5;
    xk(2) = (xk(1) - dp) + rho * g * hi * 1e-5;
    xk(3) = cv * uk(3) * ...
            (sqrt((pnom - pe) * 1e5 / rho) + alpha * (pci));
end

% Pump curve approximation
function h = pump(u, q, p1, p2, p3)
    h = p1 * u + p2 * u * q + p3 * q^2;
end
