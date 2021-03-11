%% Problema 2: Metodo BT - CS para la ecuacion de Black - Scholes
% Apartado C - Calcular la vega
% V_t = 1/2 S^2 V_ss + rS V_s -rV
% con las condiciones 
% V(S,T) = max(S-K,0) (European call option) 
% u(0, t) = 0 % Cond. Dirichlet
% u(Smax, t) = Smax - K e^{-r (T-t)}.
% en el dominio 0 <= S <= Smax, 0 <= t <= T
%% Limpiamos
clear('all')
close('all')
clc

%% Parametros del problema
r = 0.2; % Tipo de interes
Nt = 1000; % Numero de puntos en la malla temporal
Ns = 100; % Numero de puntos en la malla de precios y volatilidades
Smax = 20; % Precio maximo del subyacente
Smin = 0; % Precio minimo del subyacente
smax = 1; % Volatilidad maxima del subyacente
smin = 0; % Volatilidad minima del subyacente
T = 1; % Periodo de madurez del contrato
E = 10; % Precio de ejecución en la madurez

dt = T / Nt; % Paso temporal
ds = (Smax-Smin) / Ns; % Paso espacial

tau = linspace(0, T, Nt+1); % tau = representa el tiempo hasta la madurez
s = linspace(Smin, Smax, Ns+1);
sigma = linspace(smin, smax, Ns); % Volatilidad del subyacente
% Numeros de Courant
lambda1 =  dt / ds^2;
lambda2 = dt / ds;

%% 2. Inicializamos la matriz solucion
V = zeros(Ns+1,Nt+1,length(sigma));

%% Condiciones inicial y contorno
for i = 1:length(sigma)
    V(:,1,i) = max(s-E,0); % Condición final en European call option
    V(Ns+1,:,i) = (Smax - E * exp(-r*tau));
    V(1,:,i) = zeros(1,Nt+1);
end
%% Definicion de la matriz A
% A = I + r K
K = toeplitz([2 -1 zeros(1,Ns-3)]);
M = toeplitz([0 1 zeros(1,Ns-3)],[0 -1 zeros(1,Ns-3)]);
S = diag(s(2:end-1));
S2 = S.^2;
A = zeros(Ns-1,Ns-1,length(sigma));
for i = 1:length(sigma)
    A(:,:,i) = (1 + r*dt)*eye(Ns-1) + 1/2 * sigma(i)^2 * lambda1 * (S2*K) +...
        r/2 * lambda2 * (S * M);
end
%% Iteracion
for i = 1:length(sigma)
    for n = 1:(Nt)
        Vstar = V(2:(end-1),n,i);
        b = [(0.5*sigma(i)^2*s(2)^2*lambda1 - r/2 * s(2) * lambda2) *...
            V(1,n+1,i), zeros(1,Ns-3), (0.5*sigma(i)^2*s(Ns)^2*lambda1 +...
            r/2 * s(Ns)* lambda2) * V(Ns+1,n+1,i)]';
        Vstar = A(:,:,i) \ (Vstar +  b);
        V(:,n+1,i) = [V(1,n+1,i); Vstar; V(Ns+1,n+1,i)];
    end
end
%% Computamos vega por diferencias finitas
h = (smax-smin) / Ns; % Paso de la sigma
vega = zeros(Ns+1,Nt+1,length(sigma));
for i = 2:length(sigma)-1
    vega(:,:,i) = (V(:,:,i+1) - V(:,:,i-1) ) / h;
end
%% Dibujamos la grafica
grafica_vega = zeros(1,length(sigma));
for i = 1:length(sigma)
    grafica_vega(i) = vega(fix((Ns+1)/2),Nt+1,i);
end
figure(1);plot(sigma(2:end-1),grafica_vega(2:end-1))
xlabel('\sigma');
ylabel('Vega(\nu)');
title({'Vega para una opción europea'; 'S = S/2, T = 0, r = 0.2'});
