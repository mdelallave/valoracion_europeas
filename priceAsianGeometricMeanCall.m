function price = priceAsianGeometricMeanCall(S0,K,r,T,sigma,N)
%% priceAsianGeometricMeanCall: Price of a Asian call option on the geometric mean in the Black-Scholes model
%
%% SYNTAX:
% price = priceAsianGeometricMeanCall(S0,K,r,T,sigma,N)
%% INPUT:
% S0 : Initial value of the underlying asset
% K : Strike
% r : Risk-free interest rate
% T : Time to expiry
% sigma : Volatility
% N : Number of monitoring times
%% OUTPUT:
% price : Price of the option in the Black-Scholes model
%% EXAMPLE:
% S0 = 100; r = 0.05; K = 90; T = 2; sigma = 0.4; N = 24;
% price = priceAsianGeometricMeanCall(S0,K,r,T,sigma,N)
% M = 1e6;
% [price_MC,std_MC] = priceAsianGeometricMeanCallMC(S0,K,r,T,sigma,M,N)

r_asian = 0.5 * (r*(N+1)/N - (sigma*sigma/6) * (1 - 1/(N*N)));

sigma_asian = sigma * sqrt((2*N*N + 3*N + 1) / (6*N*N));

d_plus = (log(S0/(K*exp(-r_asian*T))) / (sigma_asian * sqrt(T))) +...
    0.5 * sigma_asian*sqrt(T);

d_minus = (log(S0/(K*exp(-r_asian*T))) / (sigma_asian * sqrt(T))) -...
    0.5 * sigma_asian*sqrt(T);

price = exp(-r*T) * (S0 * exp(r_asian*T) * normcdf(d_plus) -...
    K * normcdf(d_minus));

