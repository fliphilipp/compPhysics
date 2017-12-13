% task 3: self-consistency loop

% set up the variables
a = 0;
b = 10;
n = 3000;
h = (b-a)/n;
onevec = ones(n,1);
r = (a+h:h:b)';
Vn = 2./r; %nuclear potential: 2/r
Vh = zeros(n,1); %hartree potential: initial guess is zero
nucPot = spdiags(Vn, 0, n, n); %get to matrix form for calculation
harPot = spdiags(Vh, 0, n, n); %get to matrix form for calculation

M = spdiags([onevec -2*onevec onevec], -1:1, n, n); %matrix for poisson's equation
Mh2 = M ./ h^2;

% the loop
i = 1;
absDiffGroundE = 1;
while absDiffGroundE > 1e-5
  
  mat = -1/2 * Mh2 - nucPot + harPot; %left side of one-dimensional radial Schröinger
  [e_funct, e_values] = eig(full(mat)); %solve eigenvalue problem
  [E, iE] = min(diag(e_values)); %get smallest eigenvalue
  fr = e_funct(:,iE); %get minimal eigenvector
  rho_r = fr.^2; %probability to find electron at distance r from nucleus
  rho_r = rho_r ./ sum(h * rho_r); % normalize f(r) such that int_0^\infty f(r)^2 dr = 1
  density = rho_r ./ (4*pi*r.^2); %get density
  
  % get hartree potential
  w = -4 * pi .* r .* density .* h^2; %set up poisson's equation
  w(end) = w(end) - 1; %boundary condition
  Ur = cgs(M, w, 1e-10, 2*n); %solve 1D poisson's equation
  
  Vh = Ur ./ r; %get hartree potential
  harPot = spdiags(Vh, 0, n, n); %get to matrix form for calculation
  
  groundE_new = 2*E - sum(h * Vh .* rho_r); %integral from equation (3)
  absDiffGroundE = abs(groundE_new - groundE); %for convergence condition
  groundE = groundE_new;
  
  disp(['iteration ' num2str(i) ':  E_0 = '  num2str(groundE,5) 'a.u. = ' num2str(groundE*27.2114,5) 'eV'])
  i = i+1;
end

%% plot
Z = 2;
cfa_unscreened = Z^3 * 4 * r.^2 .* exp(-2 * Z .* r);
Z = 27/16;
cfa_opt = Z^3 * 4 * r.^2 .* exp(-2 * Z .* r);

lw = 3;
fig3 = figure(3); set(fig3, 'Position', [100, 0, 1300, 800]);
hold on
plot(r,cfa_unscreened,'color', [0.5 0.5 1],'linewidth',lw)
plot(r,cfa_opt,'color', [1 0.5 0.5],'linewidth',lw)
rho_r = rho_r ./ sum(h* rho_r);
plot(r,rho_r,'k--','linewidth',lw)

set(gca, 'fontsize',20)
ylabel('probability $\rho(r)^2$','interpreter','latex','fontsize',40)
xlabel('$r\ [a.u.]$','interpreter','latex','fontsize',40)
hartreeDescrip = ['finite difference hartree solution: $ E = ' num2str(groundE,4) '\ a.u.$'];
lgd = legend('unscreened central field: $E = -2.750\ a.u.$',...
  'optimized central field: $E = -2.848\ a.u.$',...
  hartreeDescrip);
set(lgd, 'interpreter','latex', 'fontsize', 30);
xlim([0 5])

saveas(fig3,'task3.png')

