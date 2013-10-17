%
% Mesh Generation -- SME5827 -- ICMC-USP
% Author: Afonso Paiva -- apneto@icmc.usp.br
% Date: 2013/08/30
%
% Gradient in Polar Coordinates
%

clear;
n = 32;
r = linspace(1,3,n);
theta = linspace(0,pi,n);
dr = r(2)-r(1);
dt = theta(2) - theta(1);

% Pontos no grid
[R,THETA] = meshgrid(r,theta);

% Mapeamento
X = R.*cos(THETA);
Y = R.*sin(THETA);

% Avaliando a funcao no dom. fisico
% func = sin(THETA) + exp(-0.5*R.^2); <==>
func = Y./sqrt(X.^2 + Y.^2) + exp(-0.5*X.^2-0.5*Y.^2);cd Dcd Dco

cosij = @(i,j) cos(THETA(i,j));
senij = @(i,j) sin(THETA(i,j));

% Diferenca finita centrada em r
dfcr  = @(i,j) (func(i,j+1)-func(i,j-1))/(2*dr);
% Diferenca finita centrada em theta
dfct  = @(i,j) (func(i+1,j)-func(i-1,j))/(2*dt);

% Diferencas finitas progressiva/regressiva em r
dffr  = @(i,j) (func(i,j+1)-func(i,j))/dr;
dfbr  = @(i,j) (func(i,j)-func(i,j-1))/dr;

% Diferencas finita progressiva/regressiva em theta
dfft  = @(i,j) (func(i+1,j)-func(i,j))/dt;
dfbt  = @(i,j) (func(i,j)-func(i-1,j))/dt;

fx = zeros(n);
fy = zeros(n);

for i=2:n-1
    for j=2:n-1
        fx(i,j) = cosij(i,j)*dfcr(i,j) - (senij(i,j)/R(i,j))*dfct(i,j);
        fy(i,j) = senij(i,j)*dfcr(i,j) + (cosij(i,j)/R(i,j))*dfct(i,j);
    end
end

for k=2:n-1
        fx(k,1) = cosij(k,1)*dffr(k,1) - (senij(k,1)/R(k,i))*dfft(k,1);
        fx(k,n) = cosij(k,n)*dfbr(k,n) - (senij(k,n)/R(k,n))*dfbt(k,n);
        fy(1,k) = senij(1,k)*dffr(1,k) + (cosij(1,k)/R(1,k))*dfft(1,k);
        fy(n,k) = senij(n,k)*dfbr(n,k) + (cosij(n,k)/R(n,k))*dfbt(n,k);
end

% Plotando
axis equal;
hold on
for i=1:n
    plot(X(i,:),Y(i,:),'k');
end
for j=1:n
    plot(X(:,j),Y(:,j),'k');
end
pcolor(X,Y,func);
drawnow, pause;
quiver(X,Y,fx,fy,'m');
hold off