function [price_MC_arithmetic,stdev_MC_arithmetic] = ...
    priceAsianArithmeticMeanCallMC_controlVariate(S0,K,r,T,sigma,M,N)
%% priceAsianArithmeticMeanCall_controlVariate: Control variate variance reduction for the
% price of a Asian call option on the arithmetic mean in the Black-Scholes model
%
%% SYNTAX:
% [price_MC,stdev_MC] = priceAsianArithmeticMeanCallMC_controlVariate(S0,K,r,T,sigma,M,N)
%
%% INPUT:
% S0 : Initial value of the underlying asset
% K : Strike
% r : Risk-free interest rate
% T : Time to expiry
% sigma : Volatility
% M : Number of simulations
% N : Number of observations
%
%% OUTPUT:
% price_MC : MC estimate of the price of the option in the Black-Scholes model
% stdev_MC : MC estimate of the standard deviation
%
%% EXAMPLE:
% S0 = 100; K = 90; r = 0.05; T = 2; sigma = 0.4;
% M = 1e5; N = 24;
% [price_MC,stdev_MC] = priceAsianArithmeticMeanCallMC_controlVariate(S0,K,r,T,sigma,M,N)
%

%% Generate M x N samples from N(0,1)
X = randn(M,N);
%% Simulate M trajectories in N steps
deltaT = T/N;
e = exp((r-0.5*sigma*sigma)*deltaT + sigma*sqrt(deltaT)*X);
S = cumprod([S0*ones(M,1) e],2);
%% Payoff for each trajectory in an arithmetic mean call
payoff_arithmetic = max(mean(S(:,2:end),2) - K,0);
%% Payoff for each trajectory in a geometric mean call
% payoff_GM = max(prod(S(:,2:end),2).^(1/N) - K,0); % may cause
% over/underflows
payoff_geometric = max(exp(mean(log(S(:,2:end)),2)) - K,0); % alternative
%% MC estimate of the price and the error of the arithmetic mean call
discountFactor = exp(-r * T);
price_MC_arithmetic = discountFactor * mean(payoff_arithmetic);
stdev_MC_arithmetic = discountFactor / sqrt(M) * std(payoff_arithmetic);
%% MC estimate of the price and the error of the geometric mean call
price_MC_geometric = discountFactor * mean(payoff_geometric);
stdev_MC_geometric = discountFactor * std(payoff_geometric) / sqrt(M);
%% Compute the price of the geometric mean call using the closed formula
price_geometric = priceAsianGeometricMeanCall(S0,K,r,T,sigma,N);
%% Control variate estimation
covariance = cov(payoff_arithmetic,payoff_geometric); % Var-Cov matrix
rho = corrcoef(payoff_arithmetic,payoff_geometric); % correlation matrix
price_MC = price_MC_arithmetic - covariance(1,2)/var(payoff_geometric) *...
    (price_MC_geometric - price_geometric);
stdev_MC = sqrt(1 - rho(1,2)*rho(1,2)) * stdev_MC_arithmetic; 

