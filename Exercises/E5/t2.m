a = 0;
b = 5;
n = 1000;
h = (b-a)/n;
onevec = ones(n,1);
r = (a+h:h:b)';

% second derivative matrix
mat = spdiags([onevec -2*onevec onevec], -1:1, n, n) ./h^2;

% the discretized left part of radial Schrödinger
mat = -1/2 * mat - spdiags(2./r, 0, n, n);

% get eigenvalues and minimal eigenvector
[e_funct, e_values] = eig(full(mat));
[E, iE] = min(diag(e_values));
disp(['ground state energy: ' num2str(E)])
fr = e_funct(:,iE);

% the wave function according to our ansatz
psi = 1/sqrt(4*pi) * fr ./ r;

% plot
fig2 = figure(2); set(fig2, 'Position', [100, 10, 800, 700]);
subplot(2,1,1)
plot(r,fr)
title('f(r)','fontsize',25)
subplot(2,1,2)
plot(r,psi)
title('wave function','fontsize',25)
xlabel(['$E_0 = ' num2str(E) '$'],'interpreter','latex','fontsize',20)

saveas(fig2, 'task2.png')

