function [price_MC,stdev_MC] = priceBasket2CallMC(S1_0,S2_0,c1,c2,K,r,T,sigma1,sigma2,rho,M)
%% priceBasket2CallMC: Price of a call option on a 2 asset basket in the Black-Scholes model
%
%% SYNTAX:
% [price_MC,stdev_MC] = priceBasket2CallMC(S1_0,S2_0,c1,c2,K,r,T,sigma1,sigma2,rho,M)
%
%% INPUT:
% S1_0 : Initial value of the underlying asset I
% S2_0 : Initial value of the underlying asset II
% c1 : coefficient of asset I in the basket
% c2 : coefficient of asset II in the basket
% K : Strike
% r : Risk-free interest rate
% T : Time to expiry
% sigma1 : Volatility of asset I
% sigma2 : Volatility of asset II
% rho : Correlation between the asset log-returns
% M : Number of simulations
%
%% OUTPUT:
% price_MC : MC estimate of the price of the option in the Black-Scholes model
% stdev_MC : MC estimate of the standard deviation
%
%% EXAMPLE:
% S1_0 = 100; c1 = 0.4; sigma1 = 0.2;
% S2_0 = 200; c2 = 0.3; sigma2 = 0.4;
% rho = 0.5;
% K = 90; r = 0.05; T = 2;
% M = 1e6;
% [price_MC,stdev_MC] = priceBasket2CallMC(S1_0,S2_0,c1,c2,K,r,T,sigma1,sigma2,rho,M)
% price = priceBasket2Call(S1_0,S2_0,c1,c2,K,r,T,sigma1,sigma2,rho)
%

%% Generate M x 2 independent samples from N(0,1)
X1 = randn(M,1);
X2 = randn(M,1);
%% Generate M x 2 correlated samples from N(0,rho)
Z1 = X1;
Z2 = rho*X1 + sqrt(1 - rho^2)*X2;
%% Simulate M trajectories in one step
S1_T = S1_0 * exp((r-(1/2)*sigma1^2)*T + sigma1*sqrt(T)*Z1);
S2_T = S2_0 * exp((r-(1/2)*sigma2^2)*T + sigma2*sqrt(T)*Z2);
%% Compute the payoff
basketValue = c1*S1_T + c2*S2_T;
payoff = max(basketValue - K, 0);
%% MC estimate of the price and the error of the option
discountFactor = exp(-r*T);
price_MC = discountFactor * mean(payoff);
stdev_MC = discountFactor * std(payoff) / sqrt(M);