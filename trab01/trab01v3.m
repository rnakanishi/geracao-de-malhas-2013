%
%   Geracao de Malhas - SME5827
%   Rafael Umino Nakanishi
%   
%   Equacao do calor em coordenadas esfericas
%

clear;
n = 32;
frames = 60;
theta = linspace(0, 2*pi,n);
phi = linspace(0,pi/4,n);
rho = 1;
h = 1e-5;

dp = phi(2) - phi(1);
dt = theta(2) - theta(1);

[PHI, THETA] = meshgrid(phi, theta);

% Condicoes iniciais
u = ones(n,n);
u(:,1) = 0;
u(:,n) = 10;

sinphi = @(i,j) sin(PHI(i,j)+pi/4);
cosphi = @(i,j) cos(PHI(i,j)+pi/4);
halfsinp = @(i,j) (sinphi(i,j+1)+sinphi(i,j))/2; % half sin plus
halfsinm = @(i,j) (sinphi(i,j+1)+sinphi(i,j))/2; % half sin minus

u0=u;
unew = u;
aux = u;
handle = figure;
colorbar;
for t=1:50
    for i=2:n-1
        for j=2:n-1
            % Discretizacao do espaco - Diferencas finitas centradas
            % Theta
            dfct = (1/(dt^2*sinphi(i,j)^2))*(u(i-1,j) - 2*u(i,j) + u(i+1,j));
            % Phi
            dfcp = (1/(dp^2*sinphi(i,j)))*( halfsinm(i,j)*u(i,j-1) - (halfsinm(i,j)-halfsinp(i,j))*u(i,j) + halfsinp(i,j)*u(i,j+1) );
            % Metodo de Euler
            unew(i,j) = u(i,j) + h*(dfct+dfcp);
        end
    end
    u = unew;
    [X Y Z] = sph2cart(THETA, PHI, 1);

    colormap('hot');
    surf(X,Y,Z,u);
    if mod(t,5)==0 
        saveas(handle,['time-' int2str(t)],'eps');
    end
    
    pause(0.01);
end