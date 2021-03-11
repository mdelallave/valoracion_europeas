function price = priceBasket2Call(S1_0,S2_0,c1,c2,K,r,T,sigma1,sigma2,rho)
%% priceBasket2CallMC: Price of a call option on a 2 asset basket in the Black-Scholes model
%
%% SYNTAX:
% price = priceBasket2Call(S1_0,S2_0,c1,c2,K,r,T,sigma1,sigma2,rho)
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
%
%% OUTPUT:
% price: Price of the option in the Black-Scholes model
%
%% USES:
% integrand_priceBasket2Call.m
%
%% EXAMPLE:
% S1_0 = 100; c1 = 0.4; sigma1 = 0.2;
% S2_0 = 200; c2 = 0.3; sigma2 = 0.4;
% rho = 0.5;
% K = 90; r = 0.05; T = 2;
% price = priceBasket2Call(S1_0,S2_0,c1,c2,K,r,T,sigma1,sigma2,rho)
% M = 1e6;
% [price_MC,stdev_MC] = priceBasket2CallMC(S1_0,S2_0,c1,c2,K,r,T,sigma1,sigma2,rho,M)
%

%% External integrand
externalIntegrand = @(x1)(normpdf(x1).*integrand_priceBasket2Call(x1,...
    S1_0,S2_0,c1,c2,K,r,T,sigma1,sigma2,rho));
%% Pricing formula
discount_factor = exp(-r*T);
R = 10.0;
price = discount_factor*(integral(externalIntegrand,-R,R,'RelTol',1e-6));