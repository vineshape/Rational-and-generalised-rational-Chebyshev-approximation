function [p, q, g, Err] = RationalBisection(f, n, m, T, eps)

% f - function that we approximate
% n - degree of the numerator
% m - degree of the denominator (constant term is fixed at 1, i.e. b0 = 1)
% T - domain
% eps - precision variable

% p - coefficient for the numerator
% q - coefficients for the denominator
% g - p/q
% Err - error (f - p/q)

% Vandermonde matrix for numerator
Vn = zeros(numel(T),n);
for i = 1:numel(T)
    degn = 0:n;
    Vn(i,degn+1) = (T(i)).^(degn);      %change T(i) appropriately for any other basis functions
end

% Vandermonde matrix for denominator
Vm = zeros(numel(T),m);
for j = 1:numel(T)
    degm = 1:m;
    Vm(j,degm) = (T(j)).^(degm);       %change T(i) appropriately for any other basis functions
end

% bisection method
L = 0;                        % lower bound
U = max(abs(f(T)));           % upper bound
u = 0;

% updating lower and upper limits of the bisection interval
while U-L > eps
    u = (U+L)/2;
    if checkValue(f, n, m, T, u, Vn, Vm)
        U=u;
    else
        L=u;
    end
end

% creating the model
[p, q, ~] = createModel(f, n, m, T, U, Vn, Vm);

% compute numerator and denominator seperately
G1 = Vn*p;
G2 = Vm*q;

% fraction of the two polynomials
g = G1 ./ (1+G2);

% errors or deviations
Err = zeros(numel(T),1);
F = f(T);
for r = 1:numel(T)
    Err(r) = F(r) - g(r);
end

end

% checking the feasibility
function C = checkValue(f, n, m, T, z, Vn, Vm)
    [~,~,u] = createModel(f, n, m, T, z, Vn, Vm);
    C = (u <= 0.000001);
end


% create model
function [p, q, z] = createModel(f, n, m, T, z, Vn, Vm)

% original function values
F = f(T);

% to store f+z values
F1 = zeros(numel(T),1);
for k = 1:numel(T)
    F1(k) = F(k)+z;
end

% to store f-z values
F2 = zeros(numel(T),1);
for l = 1:numel(T)
    F2(l) = F(l)-z;
end  

% LHS matrix for coefficients in the inequalities
A = [Vn -F1.*Vm -ones(numel(T),1); -Vn F2.*Vm -ones(numel(T),1); zeros(numel(T),n+1) -Vm zeros(numel(T),1); zeros(numel(T),n+1) zeros(numel(T),m) -ones(numel(T),1)];

% to store delta values
D = zeros(numel(T),1);
for i = 1:numel(T)
    D(i) = 0.9;
end

% RHS matrix for coefficients in the inequalities
b = [F1; -F2; D; ones(numel(T),1)];

% objective function
F = [zeros(1,n+m+1) 1];

options = optimoptions('linprog','Display','none');

% solving the linear programming problem
[x] = linprog(F, A, b, [], [], [], [], options);


% store the results 
p=x(1:n+1);
q=x(n+2:m+n+1);
z=x(n+m+2);

end


























