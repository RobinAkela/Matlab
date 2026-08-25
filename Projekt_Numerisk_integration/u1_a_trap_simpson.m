function u1_a_trap_simpson

%% Parametrar
R      = 3;          % radie
alpha  = 0.9;        % formfaktor
a      = 0;          % nedre integrationsgräns
b      = R;          % övre integrationsgräns

nFixed = [30 60 120];    % n-värden som ska rapporteras i tabell
nStart = 30;         % start-n för steglängdshalvering
tol    = 0.5e-4;     % mål: 4 korrekta decimaler i volymen

stepMethod = @simpson1d;      % metod i steglängdshalveringen (@trap1d / @simpson1d)

%% 1) TABELL: TRAPETS VS SIMPSON FÖR GIVNA n
fprintf("=== Uppgift 1(a) – fasta n ===\n");
for n = nFixed
    Vtrap = volume1D(@trap1d,   a, b, n, R, alpha);
    Vsimp = volume1D(@simpson1d,a, b, n, R, alpha);
    h = (b-a)/n;

    fprintf("n = %3d, h = %.4f\n", n, h);
    fprintf("  Trapets : %.15f\n", Vtrap);
    fprintf("  Simpson : %.15f\n\n", Vsimp);
end

% 2) STEGLÄNGDSHALVERING MED VALFRI METOD
fprintf("=== Steglängdshalvering (%s) ===\n", func2str(stepMethod));

volFun = @(n) volume1D(stepMethod, a, b, n, R, alpha);
[Vfinal, nFinal] = step_halving(volFun, nStart, tol, a, b);

fprintf("\nSlutlig volym (%s, 4 dec): V ≈ %.10f (n = %d)\n", ...
        func2str(stepMethod), Vfinal, nFinal);
end

%% GENERELLA HJÄLPFUNKTIONER

function V = volume1D(quadFun, a, b, n, R, alpha)
% Beräknar volymen 2*pi*∫ integrand(r) dr med vald kvadratur
I = quadFun(@(r) integrand(r,alpha,R), a, b, n);
V = 2*pi*I;
end

function [V2,n] = step_halving(volFun, n0, tol, a, b)
% Halverar h (dubblerar n) tills |V_{2h} - V_h| < tol
n  = n0;
V2 = volFun(n);

while true
    n  = 2*n;
    V1 = V2;
    V2 = volFun(n);
    h  = (b-a)/n;
    errEst = abs(V2 - V1);

    fprintf("n = %5d, h = %.6f, V ≈ %.10f, |Δ| ≈ %.3e\n", ...
            n, h, V2, errEst);

    if errEst < tol
        break;
    end
end
end

function val = integrand(r, alpha, R)
% integrand = [g(R) - g(r)] * r
val = (g_glass(R,alpha) - g_glass(r,alpha)) .* r;
end

function g = g_glass(r, alpha)
R = 3;
g = 6 .* r.^3 .* exp(-r) .* (3 + 2*sin(alpha*R)) ./ (3 + 2*sin(alpha*r));
end

function I = trap1d(f, a, b, n)
x = linspace(a,b,n+1);
h = (b-a)/n;
y = f(x);
I = h * (0.5*y(1) + sum(y(2:end-1)) + 0.5*y(end));
end

function I = simpson1d(f, a, b, n)
if mod(n,2) ~= 0
    error('Simpson kräver jämnt n.');
end
x = linspace(a,b,n+1);
h = (b-a)/n;
y = f(x);
I = h/3 * (y(1) + y(end) + 4*sum(y(2:2:end)) + 2*sum(y(3:2:end-1)));
end
