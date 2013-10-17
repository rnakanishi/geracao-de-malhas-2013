
dom = 'm-ttm';

[errj, h] = jacobi;
saveas(h, ['figures/' dom '-jacobi.eps'], 'eps');
[errg, h] = gauss;
saveas(h, ['figures/' dom '-gauss.eps'], 'eps');
% [errt, h] = thomas;
% saveas(h, ['figures/' dom '-thomas.eps'], 'eps');
[errs, h] = sor;
saveas(h, ['figures/' dom '-sor.eps'], 'eps');

h=figure;

plot(1:length(errj), errj, 'blue'); hold on;
plot(1:length(errg), errg, 'red'); hold on;
% plot(1:length(errt), errt, 'black'); hold on;
plot(1:length(errs), errs, 'green'); hold on;

legend('Jacobi', 'Gauss', 'SOR', 'Location', 'NorthEast')
xlabel('Iterations')
ylabel('Error')
% saveas(h,['figures/' dom '-error.eps'], 'epsc')