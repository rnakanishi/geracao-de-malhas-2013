%
%   Geracao de Malhas - SME5827
%   Rafael Umino Nakanishi
%   
%   Equacao do calor em coordenadas esfericas
%

clear;
nphi = 32;
nthe = nphi;
frames = 60;
theta = linspace(0, 2*pi,nthe);
phi = linspace(0.1,pi/4,nphi);
rho = 1;

h = 0.2;
dp = phi(2) - phi(1);
dt = theta(2) - theta(1);

[THETA, PHI] = meshgrid(theta, phi);
%PHI = PHI';

% Condicoes de contorno
u = ones(nphi,nthe);
u(1,:) = 10;
u(nthe,:) = 0;

%Mapeamento
X = sin(PHI)*cos(THETA);
Y = sin(PHI)*sin(THETA);
Z = cos(PHI); 

sinphi = @(i,j) sin(PHI(i,j)+pi/4);
cosphi = @(i,j) cos(PHI(i,j)+pi/4);

% Diferencas finitas centradas
dfct2 = @(i,j) (u(i-1,j) - 2*u(i,j) + u(i+1,j))/dt^2;
dfcp1 = @(i,j) (u(i,j+1) - u(i,j-1))/2*dp;
dfcp2 = @(i,j) (u(i,j-1) - 2*u(i,j) + u(i,j+1))/dp^2;

%Diferencas finitas progressiva/regressiva
dfpt2 = @(i,j) ( - 2*u(i,j) + u(i+1,j))/dt^2;
dfpp1 = @(i,j) (u(i,j+1) )/2*dp;
dfpp2 = @(i,j) ( - 2*u(i,j) + u(i,j+1))/dp^2;

dfrt2 = @(i,j) (u(i-1,j) - 2*u(i,j))/dt^2;
dfrp1 = @(i,j) ( - u(i,j-1))/2*dp;
dfrp2 = @(i,j) (u(i,j-1) - 2*u(i,j))/dp^2;


unew = u;
    for i=2:nphi-1
        for j=2:nthe-1
            unew(i,j) = (1/sinphi(i,j)^2)*dfct2(i,j) + (cosphi(i,j)/sinphi(i,j))*dfcp1(i,j) + dfcp2(i,j); 
        end 
    end
%     for i=2:nphi-1
%             unew(i,1) = (1/sinphi(i,1)^2)*dfpt2(i,1) + (cosphi(i,1)/sinphi(i,1))*dfpp1(i,1) + dfpp2(i,1);
%             unew(i,nphi) = (1/sinphi(i,nphi)^2)*dfrt2(i,nphi) + (cosphi(i,nphi)/sinphi(i,nphi))*dfrp1(i,nphi) + dfrp2(i,nphi);
%             unew(1,i) = (1/sinphi(1,i)^2)*dfpt2(1,i) + (cosphi(1,i)/sinphi(1,i))*dfpp1(1,i) + dfpp2(1,i); 
%             unew(nphi,i) = (1/sinphi(nphi,i)^2)*dfrt2(nphi,i) + (cosphi(nphi,i)/sinphi(nphi,i))*dfrp1(nphi,i) + dfrp2(nphi,i); 
%     end