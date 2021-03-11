function [price_MC, stdev_MC] = priceEuropeanOption_MC(S0,r,T,sigma,payoff,M)
%% priceEuropeanOption_MV: MC estimate of the price of a European option in the BS model
%
%% SYNTAX:
% [price_MC, stdev_MC] = priceEuropeanOption_MC(S0,r,T,sigma,payoff,M)
%
%% INPUT:
% S0 : Initial value of the underlying asset
% r : Risk-free interest rate
% T : Time to expiry
% sigma : Volatility
% payoff : Handle to the function of ST that specifies the payoff
% M : Number of simulations
%
%% OUTPUT:
% price_MC : MC estimate of the price of the option in the Black-Scholes model
% stdev_MC : MC estimate of the error in the MC price estimate
%
%% EXAMPLE 1:
% S0 = 100; K = 90; r = 0.05; T = 2; sigma = 0.4;
% payoff = @(ST)(max(ST-K,0)); % payoff of a European call option
% M = 1e6;
% [price_MC, stdev_MC] = priceEuropeanOption_MC(S0,r,T,sigma,payoff,M)
% price = priceEuropeanOption(S0,r,T,sigma,payoff) % 30.7619
%% EXAMPLE 2:
% S0 = 100; K = 90; r = 0.05; T = 2; sigma = 0.4;
% payoff = @(ST)(max(K-ST,0)); % payoff of a European put option
% M = 1e6;
% [price_MC, stdev_MC] = priceEuropeanOption_MC(S0,r,T,sigma,payoff,M)
% price = priceEuropeanOption(S0,r,T,sigma,payoff) % 12.1973

%% generate M samples from N(0,1)
X = randn(M,1);
%% simulate M trajectories in one step
ST =  S0*exp((r - 0.5*sigma^2)*T + sigma*sqrt(T)*X);
%% MC estimate
discountFactor = exp(-r*T);
price_MC = discountFactor * mean(payoff(ST));
stdev_MC = discountFactor * std(payoff(ST)) / sqrt(M);