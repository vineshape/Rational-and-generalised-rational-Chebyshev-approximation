clear;
close all;

% function that we need to approximate
f = @(T) sqrt(abs(T-0.25));

% some more functions to try
%abs(T);
%cos(2*pi*T + 0.25) + sin(2*pi*T);
%sqrt(abs(1 - T.^3));
%sqrt(abs(1 - 4*T));
%sum(sqrt(max(1-(1.2.^i).*T , 1./(1-(1.2.^i)).*(1-(1.2.^i).*T))));
%cos(T);
%gamma(T-1);
%sqrt(abs(T-0.25));
%sin(T);

n = 4;
m = 4;
eps = 0.00001;
T = -1:0.01:1;  %discretised domain

%T = -pi:(pi/2)/100:pi;

% legend
fname = '$\sqrt{\abs{T- 2.5}}$';
appname = 'rational type (4,4)';

% computing the approximation
[p, q, g, Err] = RationalBisection(f, n, m, T, eps);

% maximal deviation
k = max(Err)

% coefficients of the numerator starting with the constant
disp('coefficients of the numerator:');
disp(p);

% coefficients of the denominator, constant term is fixed at 1
disp('coefficients of the denominator:');
disp(q);


% graph of the approximation
figure;
plot(T,f(T),'linewidth',3);
hold on
plot(T,g,':','linewidth',3);
L = legend(fname, appname);
set(L,'Interpreter','latex');
set(gca,'FontSize',26);
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPositionMode','Manual');


nameit = 'Approximation';
saveas(gcf,nameit,'pdf');
print(gcf,'-dpdf','-fillpage',nameit);


%graph of the deviations
figure;
plot(T,Err,'linewidth',3);
hold on
yline(k,'--r','linewidth',3);
hold on
yline(-k,'--r','linewidth',3);
ylabel('Approximation error')
set(gca,'FontSize',26);
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperPositionMode','Manual');

nameit = 'Error curve';
saveas(gcf,nameit,'pdf');
print(gcf,'-dpdf','-fillpage',nameit);
