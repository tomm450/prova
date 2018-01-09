a = load('5.mat');  % 5 deg
b = load('10.mat'); % 10 deg

% plot differenza 
semilogy(a.n(2:end),abs(diff(a.f(:,1))),'r',b.n(2:end),abs(diff(b.f(:,1))))
legend('5 deg','10 deg');
xlabel('division along chord');
ylabel('Cl_i - Cl_{i-1}')

figure; plot(a.n,a.t,'r',b.n,b.t);
xlabel('division along chord');
ylabel('Computational time [s]')