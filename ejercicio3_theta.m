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

r = 0.2; % Tipo de interes
sigma = 0.25; % Volatilidad del subyacente
Nt = 1000; % Numero de puntos en la malla temporal
Ns = 100; % Numero de puntos en la malla de precios
Smax = 20; % Precio maximo del subyacente
Smin = 0; % Precio minimo del subyacente
T = [0.8 1.2]; % Periodo de madurez del contrato
E = 10; % Precio de ejecución en la madurez

dt = T / Nt; % Paso temporal
ds = (Smax-Smin) / Ns; % Paso espacial
% tau = representa el tiempo hasta la madurez
for i = 1:length(T)
    tau(i,:) = linspace(0, T(i), Nt+1);
end

s = linspace(Smin, Smax, Ns+1);
% Numeros de Courant
lambda1 =  dt / ds^2;
lambda2 = dt / ds;

%% Inicializamos la matriz solucion
V = zeros(Ns+1,Nt+1,length(T));

%% Condiciones inicial y contorno
for i = 1:length(T)
    V(:,1,i) = max(s-E,0); % Condición final en European call option
    V(Ns+1,:,i) = (Smax - E * exp(-r*tau(i,:)));
    V(1,:,i) = zeros(1,Nt+1);
end
%% Definicion de la matriz A
% A = I + r K
K = toeplitz([2 -1 zeros(1,Ns-3)]);
M = toeplitz([0 1 zeros(1,Ns-3)],[0 -1 zeros(1,Ns-3)]);
S = diag(s(2:end-1));
S2 = S.^2;
A = zeros(Ns-1,Ns-1,length(T));
for i = 1:length(T)
    A(:,:,i) = (1 + r*dt(i)) * eye(Ns-1) + ...
        1/2 * sigma^2 * lambda1(i) * (S2*K) + r/2 * lambda2(i) * (S * M);
end
%% Iteracion
for i = 1:length(T)
    for n = 1:(Nt)
        Vstar = V(2:(end-1),n,i);
        b = [(0.5*sigma^2*s(2)^2*lambda1(i) - r/2 * s(2) * lambda2(i)) *...
            V(1,n+1,i), zeros(1,Ns-3), (0.5*sigma^2*s(Ns)^2*lambda1(i) +...
            r/2 * s(Ns)* lambda2(i)) * V(Ns+1,n+1,i)]';
        Vstar = A(:,:,i) \ (Vstar +  b);
        V(:,n+1,i) = [V(1,n+1,i); Vstar; V(Ns+1,n+1,i)];
    end
end
%% Theta
h = T(2) - T(1); % Paso de la T
theta = (V(:,:,2) - V(:,:,1) ) / h;

%% Dibujamos la grafica
plot(s,theta(:,1),...
    s,theta(:,fix((Nt+1)/6)),...
    s,theta(:,fix((Nt+1)/2)),...
    s,theta(:,fix((Nt+1)/3)),...
    s,theta(:,fix((Nt+1)/1.5)),...
    s,theta(:,fix((Nt+1)/1.2)),...
    s,theta(:,Nt+1))
legend('\tau = 0','\tau = 16% T','\tau = 33% T','\tau = 50% T',...
    '\tau = 66% T','\tau = 83% T','\tau = T','Location','northwest')
xlabel('Valor del subyacente (S)');
ylabel('Theta(\Theta)');
title({'Theta para una opción europea con T = 1 ';'E = 10, \sigma = 0.25, r = 0.2'});