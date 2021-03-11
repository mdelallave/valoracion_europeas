function externalIntegrand = integrand_priceBasket2Call(x1,S1_0,S2_0,c1,c2,K,r,T,sigma1,sigma2,rho)
%% integrand_priceBasket2Call: auxiliary function for priceBasket2Call
%
%% SYNTAX:
% externalIntegrand = integrand_priceBasket2Call(x1,S1_0,S2_0,c1,c2,K,r,T,sigma1,sigma2,rho)
%
%% INPUT:
% x1 : integration variable
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
% externalIntegrand : External integrand to price of the option in the Black-Scholes model
%
%% Initialize output
externalIntegrand = zeros(size(x1));
%% Vectorize
for i = 1:length(x1)
    internalIntegrand = @(x2) internal_integrand_priceBasket2Call(x2,...
        x1(i),S1_0,S2_0,c1,c2,K,r,T,sigma1,sigma2,rho); 
    R = 10.0; % Proxy for \infty for N(0,1)
    externalIntegrand(i) = integral(internalIntegrand,-R,R,'RelTol',1.0e-6);
end