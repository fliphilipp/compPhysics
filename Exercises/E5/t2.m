a = 0;
b = 7;
n = 10000;
h = (b-a)/n;
onevec = ones(n,1);
r = (a+h:h:b)';
lwd = 3;

% second derivative matrix
mat = spdiags([onevec -2*onevec onevec], -1:1, n, n) ./h^2;

% the discretized left part of radial Schrödinger
nucPot = 1./r;
mat = -1/2 * mat - spdiags(nucPot, 0, n, n);

% get eigenvalues and minimal eigenvector
% min e_value is ground state energy, minimal eigenvector
% min e_funct is non-normalized f(r)
[e_funct, e_values] = eig(full(mat));
[E, iE] = min(diag(e_values));
disp(['ground state energy: ' num2str(E)])
fr = e_funct(:,iE);

% normalize f(r) such that int_0^\infty f(r)^2 dr = 1
fr2 = fr.^2;
fr2 = fr2 ./ sum(fr2*h);
fr = sqrt(fr2);

% the wave function according to our ansatz
psi = 1/sqrt(4*pi) * fr ./ r; % the equation from the assignment
psi2 = psi.^2;

% check cusp condition
dfdr0 = (fr(2) - fr(1)) / h;
psi02 = 1 / (4*pi) * dfdr0^2;

% plot f(r)^2
fig2 = figure(2); set(fig2, 'Position', [100, 0, 800, 800]);
subplot(2,1,1)
hold on
plot(r,fr2,'linewidth',lwd)
set(gca, 'fontsize',16)
ylabel('$f(r)^2$','interpreter','latex','fontsize',35)
xlabel('$r\ [a.u.]$','interpreter','latex','fontsize',35)

% plot psi^2
subplot(2,1,2)
plot(r,psi2,'linewidth',lwd)
set(gca, 'fontsize',16)
xlabel(['$E_0 = ' num2str(E) ' H_a = ' num2str(E*27.21138602) 'eV$'],'interpreter','latex','fontsize',35)
hold on
psi_theo = 1/sqrt(pi) .* exp(-r);
psi2_theo = psi_theo.^2;

plot(r,psi2_theo,'--','linewidth',lwd)
scatter(0,psi02,80,'black','filled')
set(gca, 'fontsize',16)
lgd = legend('numeric result for $|\psi(r)|^2$', 'analytical result for $|\psi(r)|^2$',...
  'cusp condition: $\lim_{r\rightarrow 0}\left(|\psi(r)|^2\right) = \frac{1}{4\pi} \left(\frac{\partial f(0)}{\partial r}\right)^2$');
set(lgd, 'interpreter','latex', 'fontsize', 20);
Edescrip = ['$E_0 = ' num2str(E) ' H_a = ' num2str(E*27.21138602) 'eV$'];
text(2.5,0.2,Edescrip,'interpreter','latex','fontsize',20)
ylabel('$ | \psi |^2$','interpreter','latex','fontsize',35)
xlabel('$r\ [a.u.]$','interpreter','latex','fontsize',35)

saveas(fig2, 'task2.png')
