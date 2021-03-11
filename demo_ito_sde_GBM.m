function demo_ito_sde_GBM
%% demo_ito_sde_GBM: Itô  stochastic differential equation for GBM
%
%   Simulate Itô process S(t) in [t0, t0+T]
% 
%   SDE:    dS(t) = a(S(t), t)*dt + b(S(t), t)*dW(t)
%
%                   a(S(t), t) = mu*S(t):    Drift term
%                   b(S(t), t) = sigma*S(t): Difussion
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
S0 = 100.0; % Initial value of the process

%% SDE for arithmetic Brownian Motion: BM(mu,sigma)
mu = 0.1;
a = @(S_t,t) mu*S_t; % drift

sigma = 0.4;
b = @(S_t,t) sigma*S_t; % difussion

M = 5000; % Number of simulations

X = randn(M,N); % M Gaussian white noise trajectories


S(:,1) = S0*ones(M,1);
tic
for n = 1:N
   S(:,n+1) = S(:,n) + a(S(:,n), t(n))*dT + ...          % Predictable
                       b(S(:,n), t(n))*sqrt(dT).*X(:,n); % Innovation 
end
toc


figure(2);
nBins = 30;
ST = S(:,end);
m = median(log(ST/S0));
s = std(log(ST/S0));

histogram(ST/S0, nBins,'Normalization','pdf');

hold on
nPlot = 1000;
alpha = 4.0;
xPlot = linspace(exp(m - alpha*s), exp(m + alpha*s), nPlot);
plot(xPlot, lognpdf(xPlot, m, s),'linewidth',2.5);
hold off

figure(1);
plot(t,log(S(1:min(50,M),:)'))
hold on
plot(t, mean(log(S)),'k','linewidth',2)
plot(t, mean(log(S)) - 2*std(log(S)),'b', 'linewidth',2)
plot(t, mean(log(S)) + 2*std(log(S)),'r', 'linewidth',2)
hold off
