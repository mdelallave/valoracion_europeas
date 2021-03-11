function price = priceEuropeanCall(S0,K,r,T,sigma)
%% priceEuropeanCall: Price of a European call option in the Black-Scholes model
%
%% SYNTAX:
% price = priceEuropeanCall(S0,K,r,T,sigma)
%
%% INPUT:
% S0 : Initial value of the underlying asset
% K : Strike
% r : Risk-free interest rate
% T : Time to expiry
% sigma : Volatility
%
%% OUTPUT:
% price : Price of the option in the Black-Scholes model
%
%% EXAMPLE:
% S0 = 100; r = 0.05; K = 90; T = 2; sigma = 0.4;
% price = priceEuropeanCall(S0,K,r,T,sigma)
%
%% Value of assets at maturity
ST = @(x) S0*exp((r - 0.5*sigma*sigma)*T + sigma*sqrt(T)*x);

%% Payoff for EU call option
payoff = @(x) max(ST(x) - K, 0); % Function of x through ST(x)

%% Integrand for expected value of payoff
integrand = @(x) normpdf(x).*payoff(x);

%%
discount_factor = exp(-r*T);
R = 10; % Proxy for \infty for N(0,1) 

% Price = expected value of payoff discounted (discounted to account for
% time-value of euros)
price = discount_factor*(integral(integrand,-R,R,'RelTol',1.0e-6));
