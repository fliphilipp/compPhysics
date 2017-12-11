
%% 1d poisson according to 
%  http://farside.ph.utexas.edu/teaching/329/lectures/node62.html

clear all
close all

a = 0;
b = 7;
n = 1000;
h = (b - a) / n;
ul = 0;
uh = 1;
r = a:h:b-h;
r = [r'; b];
onevec = ones(n+1,1);

rho = 4 * r.^2 .* exp(-2.*r); % hydrogen ground state electron density
subplot(2,1,1)
plot(r,rho)
lgd = legend('$n_s(r)$ electron density');
set(lgd, 'interpreter','latex', 'fontsize', 14);

% discretize second derivative
M = spdiags([onevec -2*onevec onevec], -1:1, n+1, n+1);
w = (-4 * pi * r .* rho) .* h^2;
w(1) = w(1) - ul;
w(end) = w(end) - uh;

Ur = cgs(M, w, 1e-7, 5000);
Vh = Ur ./ r;
subplot(2,1,2)
hold on
plot(r,Ur,'g-')
plot(r,Vh)

Vh_theory = 1 ./ r - (1 + 1 ./ r) .* exp(-2*r);
plot(r, Vh_theory)
lgd = legend('$U(r)$', '$V_h$','$V_h$ (theoretical)');
set(lgd, 'interpreter','latex', 'fontsize', 14);
