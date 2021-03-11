function [vega_MC,stdev_MC] = vegaEuropeanCallMC(S0,K,r,T,h,sigma,M)
%% vegaEuropeanCallMC: Vega of a European call option in the Black-Scholes model
%
%% SYNTAX:
% [vega_MC,stdev_MC] = vegaEuropeanCallMC(S0,K,r,T,h,sigma,M)
%
%% INPUT:
% S0 : Initial value of the underlying asset
% K : Strike
% r : Risk-free interest rate
% T : Time to expiry
% h : error
% sigma : Volatility
% M : Number of simulations
%
%
%% OUTPUT:
% vega_MC : MC estimate of the vega of the option in the Black-Scholes model
% stdev_MC : MC estimate of the standard deviation
%
%% EXAMPLE:
% S0 = 100; K = 90; r = 0.05; T = 2; h = 1.0e-5; sigma = 0.4;
% M = 1e6;
% [vega_MC,stdev_MC] = vegaEuropeanCallMC(S0,K,r,T,h,sigma,M)
% vega = vegaEuropeanCall(S0,K,r,T,sigma)
%% generate M samples from N(0,1)
X = randn(M,1);
%% simulate M minus / plus trajectories in one step
sigmap = sigma*(1 + h);
sigmam = sigma*(1 - h);
STp =  S0*exp((r - 0.5* sigmap*sigmap)*T + sigmap*sqrt(T)*X);
STm =  S0*exp((r - 0.5* sigmam*sigmam)*T + sigmam*sqrt(T)*X);
%% define minus / plus payoffs
payoffp = max(STp-K,0);
payoffm = max(STm-K,0);
%% compute vega along each trajectory
vega = (payoffp - payoffm) ./ (2*sigma*h);
%% MC estimate
discountFactor = exp(-r*T);
vega_MC = discountFactor * mean(vega);
stdev_MC = discountFactor * std(vega) / sqrt(M);