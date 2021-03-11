function demo_simulacion_browniano_geometrico

S0 = 100.0;
r = 0.5;
sigma = 0.4;

T = 2; % years for simulation
N = 8; % quarterly
dT = T/N; % step size

M = 5; % number of simulations


X = randn(M,N); % X ~ N(0,1)

e = exp((r-0.5*sigma*sigma)*dT + sigma*sqrt(dT)*X); % e ~ LN
St = cumprod([S0*ones(M,1) e], 2) % simulations
figure(1); plot(0:dT:T, ST')

ST = St(:,end); % final value for each simulation

payoff_EU_call = max (ST-K,0); % EU call plain vanilla

payoff_call_arithmetic_mean_MC = ...;

discount_factor = exp(-r*T);

price_EU_call_MC = discount_factor * mean(payoff_EU_call)
stdev_EU_call_MC = discount_factor * std(payoff_EU_call)/sqrt(M)