function demo_ito_sde
%% demo_ito_sde: Itô  stochastic differential equation
%
%   Simulate Itô process Z(t) in [t0, t0+T]
% 
%   SDE:    dZ(t) = a(Z(t), t)*dt + b(Z(t), t)*dW(t)
%
%                   a(Z(t), t): Drift term
%                   b(Z(t), t): Difussion
%                        dW(t): Wiener process (standard Brownian Motion)
%

%% Simulate between [t0, t0+T]
t0 = 1.0;
T = 3.0; % Length of simulation

%% Number of steps for simulation
days_in_year = 252;
N = T*days_in_year; 
dT = T/N; % Should be as small as possible
t = t0:dT:(t0+T); % Monitoring times

%% SDE for the process
Z0 = 0.0; % Initial value of the process

%% SDE for arithmetic Brownian Motion: BM(mu,sigma)
mu = 3.0;
a = @(Z_t,t) mu; % drift

sigma = 2.0;
b = @(Z_t,t) sigma; % difussion

M = 5000; % Number of simulations

X = randn(M,N); % M Gaussian white noise trajectories


Z(:,1) = Z0*ones(M,1);
tic
for n = 1:N
   Z(:,n+1) = Z(:,n) + a(Z(:,n), t(n))*dT + ...         % Predictable
                       b(Z(:,n), t(n))*sqrt(dT)*X(:,n); % Innovation 
end
toc


figure(2);
nBins = 30;
ZT = Z(:,end);
m = mean(ZT); % mu*T
s = std(ZT);  % sigma*sqrt(T)

histogram(ZT, nBins,'Normalization','pdf');

hold on
nPlot = 1000;
alpha = 4.0;
xPlot = linspace(m - alpha*s, m + alpha*s, nPlot);
plot(xPlot, normpdf(xPlot, m, s),'linewidth',2.5);
hold off

figure(1);
plot(t,Z(1:min(50,M),:)')
hold on
plot(t, mean(Z),'k','linewidth',2)
plot(t, mean(Z) - 2*std(Z),'b', 'linewidth',2)
plot(t, mean(Z) + 2*std(Z),'r', 'linewidth',2)
hold off
