
%% 1d poisson according to 
%  http://farside.ph.utexas.edu/teaching/329/lectures/node62.html

clear all
close all

a = 0;
b = 5;
n = 1000;
h = (b - a) / n;
ul = 0;
uh = 1;
r = a+h:h:b;
r = [r'; b];
onevec = ones(n+1,1);
lwd = 3;

fig1 = figure(1); set(fig1, 'Position', [100, 10, 800, 700]);

a0 = 1;
rho_r = 4 / a0^3 * r.^2 .* exp(-2/a0*r); % hydrogen ground state electron density
rho = a0^3 / pi .* exp(-2/a0*r); % leave the 4 * pi * r.^2 term out for whatever reason
subplot(2,1,1)
plot(r,rho_r,'linewidth',lwd)
sum(rho .* h)
set(gca,'fontsize',15)
lgd = legend('$n_s(r)$ electron density');
xlabel('distance $r\ [a.u.]$','fontsize', 25, 'interpreter','latex')
ylabel('density','fontsize', 25, 'interpreter','latex')
set(lgd, 'interpreter','latex', 'fontsize', 18);

% discretize second derivative
M = spdiags([onevec -2*onevec onevec], -1:1, n+1, n+1);
w = (-4 * pi .* r .* rho);
w =  w .* h^2;
w(1) = w(1) - ul;
w(end) = w(end) - uh;

Ur = cgs(M, w, 1e-10, 2*n);
Vh = Ur ./ r;
subplot(2,1,2)
hold on
plot(r,Ur,'linewidth',lwd)
plot(r,Vh,'linewidth',lwd)
xlabel('distance $r\ [a.u.]$','fontsize', 25, 'interpreter','latex')
ylabel('potential','fontsize', 25, 'interpreter','latex')

Vh_theory = 1 ./ r - (1 + 1 ./ r) .* exp(-2*r);
Ur_theory = Vh_theory.*r;
plot(r, Ur_theory,'--','linewidth',lwd)
plot(r, Vh_theory,'--','linewidth',lwd)
lgd = legend('$U(r)$', '$V_h(r)$','$U(r)$ (theoretical)','$V_h(r)$ (theoretical)');
set(lgd, 'interpreter','latex', 'fontsize', 18, 'location', 'east');
ylim([0 1.1])

saveas(fig1, 'task1.png')

figure(123)
secdiff = nan(n+1,1);
for i = 2:n
  secdiff(i) = (Ur_theory(i-1) + Ur_theory(i+1) - 2*Ur_theory(i)) ./ (h.^2);
end
hold on
plot(r, secdiff,'linewidth',lwd)
plot(r, [nan; w(2:end-1)./h.^2; nan],'--','linewidth',lwd)


