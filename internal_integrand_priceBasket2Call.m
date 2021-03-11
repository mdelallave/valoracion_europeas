function internalIntegrand = internal_integrand_priceBasket2Call(x2,x1,S1_0,S2_0,c1,c2,K,r,T,sigma1,sigma2,rho)
% COMENTAR

S1_T = S1_0 * exp((r-(1/2)*sigma1^2)*T + sigma1*sqrt(T)*x1);
S2_T = S2_0 * exp((r-(1/2)*sigma2^2)*T + ...
    sigma2*sqrt(T)*(rho*x1 + sqrt(1 - rho^2)*x2));

payoff = max(c1*S1_T + c2*S2_T - K, 0);
internalIntegrand = normpdf(x2) .* payoff;