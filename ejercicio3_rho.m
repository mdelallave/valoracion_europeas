%% Problema 3: Metodo BT - CS para la ecuacion de Black - Scholes
% Calculo de las griegas (delta)
% V_t = 1/2 S^2 V_ss + rS V_s -rV
% con las condiciones 
% V(S,T) = max(S-K,0) (European call option) 
% u(0, t) = 0 % Cond. Dirichlet
% u(Smax, t) = Smax - K e^{-r (T-t)}.
% en el dominio 0 <= S <= Smax, 0 <= t <= T
%% Limpiamos
clearvars
close('all')
clc

%% Parametros del problema

r = [0.16 0.24]; % Tipo de interes
sigma = 0.25; % Volatilidad del subyacente
Nt = 1000; % Numero de puntos en la malla temporal
Ns = 100; % Numero de puntos en la malla de precios
Smax = 20; % Precio maximo del subyacente
Smin = 0; % Precio minimo del subyacente
T = 1; % Periodo de madurez del contrato
E = 10; % Precio de ejecución en la madurez

dt = T / Nt; % Paso temporal
ds = (Smax-Smin) / Ns; % Paso espacial

tau = linspace(0, T, Nt+1); % tau = representa el tiempo hasta la madurez
s = linspace(Smin, Smax, Ns+1);
% Numeros de Courant
lambda1 =  dt / ds^2;
lambda2 = dt / ds;

%% Inicializamos la matriz solucion
V = zeros(Ns+1,Nt+1,length(r));

%% Condiciones inicial y contorno
for i = 1:length(r)
    V(:,1,i) = max(s-E,0); % Condición final en European call option
    V(Ns+1,:,i) = (Smax - E * exp(-r(i)*tau));
    V(1,:,i) = zeros(1,Nt+1);
end
%% Definicion de la matriz A
% A = I + r K
K = toeplitz([2 -1 zeros(1,Ns-3)]);
M = toeplitz([0 1 zeros(1,Ns-3)],[0 -1 zeros(1,Ns-3)]);
S = diag(s(2:end-1));
S2 = S.^2;
A = zeros(Ns-1,Ns-1,length(r));
for i = 1:length(r)
    A(:,:,i) = (1 + r(i)*dt) * eye(Ns-1) +...
        1/2 * sigma^2 * lambda1 * (S2*K) + r(i)/2 * lambda2 * (S * M);
end
%% Iteracion
for i = 1:length(r)
    for n = 1:(Nt)
        Vstar = V(2:(end-1),n,i);
        b = [(0.5*sigma^2*s(2)^2*lambda1 - r(i)/2 * s(2) * lambda2) *...
            V(1,n+1,i), zeros(1,Ns-3), (0.5*sigma^2*s(Ns)^2*lambda1 +...
            r(i)/2 * s(Ns)* lambda2) * V(Ns+1,n+1,i)]';
        Vstar = A(:,:,i) \ (Vstar +  b);
        V(:,n+1,i) = [V(1,n+1,i); Vstar; V(Ns+1,n+1,i)];
    end
end
%% Vega
h = r(2) - r(1); % Paso de la r
rho = (V(:,:,2) - V(:,:,1) ) / h;

%% Dibujamos la grafica
plot(s,rho(:,1),...
    s,rho(:,fix((Nt+1)/6)),...
    s,rho(:,fix((Nt+1)/2)),...
    s,rho(:,fix((Nt+1)/3)),...
    s,rho(:,fix((Nt+1)/1.5)),...
    s,rho(:,fix((Nt+1)/1.2)),...
    s,rho(:,Nt+1))
legend('\tau = 0','\tau = 16% T','\tau = 33% T','\tau = 50% T',...
    '\tau = 66% T','\tau = 83% T','\tau = T','Location','northwest')
xlabel('Valor del subyacente (S)');
ylabel('Rho(\rho)');
title({'Rho para una opción europea con r = 0.2';'E = 10, T = 1, \sigma = 0.25'});