%% Volt-var Control
% Output optimal DER reactive power setpoints

clear all, close all, clc

run('vvc_setup_IEEE123s.m')

% Voltage magnitude limits (p.u.)
vmax = 1.05;
vmin = 0.95;

Y = R*-P + X*-Q + ones(N,T);
Qg = zeros(N,T);

for n = 1:N 
    Srate(n) = 2.25*max(Pg(n,:)); 
end

Qlim = zeros(N,T);
for n = 1:N
    Qlim(n,:) = sqrt(Srate(n)^2 - Pg(n,:).^2);
end

for t = 1:T
    % solving for optimal q_g using quadprog
    
    % objective variables
    f = -2*R*Q(:,t);
    H = 2*R;
    
    % voltage constraint and reactive power limits
    Eye = eye(N);
    A = [X; -X; Eye; -Eye];
    b = [-Y(:,t) + vmax*ones(N,1); Y(:,t) - vmin*ones(N,1); Qlim(:,t); Qlim(:,t)];
       
    % solve for q_g
    q_g = quadprog(H,f,A,b);
    Qg(:,t) = q_g;
end

