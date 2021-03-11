function [price_MC,stdev_MC] = priceEuropeanCallMC(S0,K,r,T,sigma,M)
%% priceEuropeanCallMC: Price of a European call option in the Black-Scholes model
%
%% SYNTAX:
% [price_MC,stdev_MC] = priceEuropeanCallMC(S0,K,r,T,sigma,M)
%
%% INPUT:
% S0 : Initial value of the underlying asset
% K : Strike
% r : Risk-free interest rate
% T : Time to expiry
% sigma : Volatility
% M : Number of simulations
%
%
%% OUTPUT:
% price_MC : MC estimate of the price of the option in the Black-Scholes model
% stdev_MC : MC estimate of the standard deviation
%
%% EXAMPLE:
% S0 = 100; K = 90; r = 0.05; T = 2; sigma = 0.4;
% M = 1e4;
% [price_MC,stdev_MC] = priceEuropeanCallMC(S0,K,r,T,sigma,M)
% M = 100*M; % lowers the MC error by a facor of 10
% [price_MC,stdev_MC] = priceEuropeanCallMC(S0,K,r,T,sigma,M)
%
%% generate M samples from N(0,1)
X = randn(M,1);
%% simulate M trajectories in one step
ST =  S0*exp((r - 0.5*sigma^2)*T + sigma*sqrt(T)*X);
%% define payoff
payoff = @(ST)(max(ST-K,0));
%% MC estimate
discountFactor = exp(-r*T);
price_MC = discountFactor * mean(payoff(ST));
stdev_MC = discountFactor * std(payoff(ST)) / sqrt(M);