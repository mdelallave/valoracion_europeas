function demo_multivariate_gaussian(rho,nPlot,nFig)
% demo_multivariate_gaussian: Multivariate Gaussian distribution
%   
%   EXAMPLE 1: 
%       rho = 0.0; % independent
%       demo_multivariate_gaussian(rho)
%       title(sprinf('Independent 2D Gaussian %s = %.2f','\rho','rho'));
%
%   EXAMPLE 2: 
%       nPlot = 50;
%       nFig = 3;
%       rho = 0.8; % positive correlated
%       demo_multivariate_gaussian(rho,nPlot,nFig)
%       title(sprinf('Independent 2D Gaussian %s = %.2f','\rho','rho'));
%
%   EXAMPLE 3: 
%       nPlot = 50;
%       nFig = 5;
%       rho = -0.8; % anticorrelated
%       demo_multivariate_gaussian(rho,nPlot,nFig)
%       title(sprinf('Independent 2D Gaussian %s = %.2f','\rho','rho'));


if (nargin < 2) 
    % default values
    nPlot = 30; 
end
if (nargin < 3)
    nFig = 1;
end

%% From 2D grid
alpha = 4;
xPlot = linspace(-alpha,alpha,nPlot);
yPlot = linspace(-alpha,alpha,nPlot);
[xPlot,yPlot] = meshgrid(xPlot,yPlot);


%% Configure pdf on grid points
Mu = zeros(2,1); % 2D Gaussian centered at the origin
Sigma = [1.0 rho;...
          rho 1.0]; % Covariance matrix

zPlot = normpdf_2D(xPlot,yPlot,Mu,Sigma);

%% Graph in 3D
figure(nFig); mesh(xPlot,yPlot,zPlot)
% figure(nFig+1); surf(xPlot,yPlot,zPlot)