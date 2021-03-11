function vega = vegaEuropeanCall(S0,K,r,T,sigma)
% vegaEuropeanCall: Vega of a European call option
%
% SYNTAX:
% vega = vegaEuropeanCall(S0,K,r,T,sigma)
%
% INPUT:
% S0 : Initial value of the underlying asset
% K : Strike
% r : Risk-free interest rate
% T : Time to expiry
% sigma : Volatility
%
% OUTPUT:
% vega : Vega of the option in the Black-Scholes model
%
% EXAMPLE:
% S0 = 100; K = 90; r = 0.05; T = 2; sigma = 0.4;
% vega = vegaEuropeanCall(S0,K,r,T,sigma)

d1 = ( log(S0./K) + (r + sigma.*sigma./2).*T ) ./ (sigma.*sqrt(T));
vega = S0.*sqrt(T).*normpdf(d1);