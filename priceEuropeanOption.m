function price = priceEuropeanOption(S0,r,T,sigma,payoff)
%% priceEuropeanOption: Price of a European option in the Black-Scholes model
%
%% SYNTAX:
% price = priceEuropeanOption(S0,r,T,sigma,payoff)
%
%% INPUT:
% S0 : Initial value of the underlying asset
% r : Risk-free interest rate
% T : Time to expiry
% sigma : Volatility
% payoff : Handle to the function of ST that specifies the payoff
%
%% OUTPUT:
% price : Price of the option in the Black-Scholes model
%
%% EXAMPLE 1:
% S0 = 100; K = 90; r = 0.05; T = 2; sigma = 0.4;
% payoff_0 = @(ST)ST; % payoff of underlying
% price = priceEuropeanOption(S0,r,T,sigma,payoff_0)
%
%% EXAMPLE 2:
% S0 = 100; K = 90; r = 0.05; T = 2; sigma = 0.4;
% payoff_1 = @(ST)(max(ST-K,0)); % payoff of a European call option
% price_EU_call = priceEuropeanOption(S0,r,T,sigma,payoff_1)
%
%% EXAMPLE 3:
% S0 = 100; K = 90; r = 0.05; T = 2; sigma = 0.4;
% payoff_2 = @(ST)(max(K-ST,0)); % payoff of a European put option
% price_EU_put = priceEuropeanOption(S0,r,T,sigma,payoff_2)
%
%       % put-call parity:
%         call = priceEuropeanOption(S0,r,T,sigma,payoff_1) + exp(-r*T)*K;
%         put = priceEuropeanOption(S0,r,T,sigma,payoff_2) + S0;
%         parity = call - put % Error about 1e-7
%
%% Value of underlying at maturity
ST = @(X) S0*exp((r - 0.5*(sigma^2))*T + sigma*sqrt(T)*X);

%% Price as expected value of payoff in risk-netrual world

discount_factor = exp(-r*T);
integrand = @(X)normpdf(X) .* payoff(ST(X));
R = 10; % Numerical support of N(0,1): [-R,R]
price = discount_factor * integral(integrand,-R,R,'RelTol',1.0e-6);