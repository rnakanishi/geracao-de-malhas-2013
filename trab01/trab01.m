%
%   Geracao de Malhas - SME5827
%   Rafael Umino Nakanishi
%   
%   Equacao do calor em coordenadas esfericas
%



clear;
n = 50;
frames = 60;
theta = linspace(0, 2*pi,n);
phi = linspace(0.1,pi/4,n);
rho = 1;

dp = phi(2) - phi(1);
dt = theta(2) - theta(1);
h = dp;

[THETA, PHI] = meshgrid(theta, phi);

% Condicoes de contorno
u = ones(n,1);
u(1) = 10;
u(n) = 0;

sinphi = @(i,j) sin(PHI(i,j)+pi/4);
cosphi = @(i,j) cos(PHI(i,j)+pi/4);

F = u;
% Iteracoes para o numero de frames
for t = 1:100
    % Sistema linear para phi
    bphi = @(i,j) (h*cosphi(i,j)/sinphi(i,j))/2;

    A = zeros(n,n);
    A(1,1) = -2;
    A(1,2) = 1 + bphi(1,2);
    A(n,n-1) = 1 - bphi(n,n-1);
    A(n,n) = -2;
    for i=2:n-1
            A(i,i-1) = 1 - bphi(i,i-1);
            A(i,i) = -2;
            A(i,i+1) = 1 + bphi(i,i+1);
    end
    A = A./h^2;
    
    F(1) = u(1) - (1/h^2 - (cosphi(1,1)/sinphi(1,1))/(2*h))*10;
    for i=2:n-1
        F(i) = u(i);
    end
    F(n) = u(n) - (1/h^2 + (cosphi(n,n)/sinphi(n,n))/(2*h))*0;
    
    u = F\A;
    u(1) = 10;
    u(50) = 0;
    
% Imprimindo imagem
[X Y Z] = sph2cart(THETA, PHI, 1);

colormap('hot');
surf(X,Y,u);
end
% Sistema linear para theta
