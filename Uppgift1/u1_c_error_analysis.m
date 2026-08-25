function u1_c_error_analysis
% Konvergensstudie för 1D- och 2D-metoder, med n och h beräknade systematiskt

%% Parametrar
R     = 3;
alpha = 0.9;
L     = 3*sqrt(2);      % storlek på kvadraten i 2D

% referens-n ("nästan exakta" lösningar)
nRef1D = 2^14;          % referens för 1D
nRef2D = 2^9;           % referens för 2D

% start-n och antal förfiningssteg
n1D_start   = 2^4;      % första n i 1D (motsvarar 16)
numLevels1D = 6;        % antal förfiningar => 16,32,...,16*2^(6-1)=512

n2D_start   = 2^3;      % första n i 2D (8)
numLevels2D = 5;        % => 8,16,32,64,128

% vilka 1D-metoder som ska jämföras
quad1D      = {@trap1d, @simpson1d};
quadNames1D = {'Trapets 1D','Simpson 1D'};
markers1D   = {'o-','s-'};

%% 1) REFERENSLÖSNINGAR
% 1D: använd Simpson med väldigt fin grid
Vref1D = volume1D(@simpson1d, R, alpha, nRef1D);

% 2D: trapets med fin grid
Vref2D = volume2D(alpha, R, L, nRef2D);

%% 2) FEL FÖR 1D-METODER (systematisk förfining)
nList1D      = zeros(1,numLevels1D);
h1           = zeros(1,numLevels1D);
err1D_all    = zeros(numel(quad1D), numLevels1D);

n = n1D_start;
for k = 1:numLevels1D
    nList1D(k) = n;
    h1(k)      = R / n;

    % beräkna fel för varje 1D-kvadratur vid detta n
    for m = 1:numel(quad1D)
        V = volume1D(quad1D{m}, R, alpha, n);
        err1D_all(m,k) = abs(V - Vref1D);
    end

    n = 2*n;  % halvera steglängden (dubblera n)
end

%% 3) FEL FÖR 2D-TRAPETS (systematisk förfining)
nList2D = zeros(1,numLevels2D);
h2      = zeros(1,numLevels2D);
err2D   = zeros(1,numLevels2D);

n = n2D_start;
for k = 1:numLevels2D
    nList2D(k) = n;
    h2(k)      = L / n;

    V = volume2D(alpha, R, L, n);
    err2D(k) = abs(V - Vref2D);

    n = 2*n;
end

%% 3b) UTSKRIFT AV NOGGRANNHETSORDNING

fprintf('\n--------------------------------------\n');
fprintf('   Beräknad noggrannhetsordning\n');
fprintf('--------------------------------------\n');

% Ordning för 1D-metoder (log-log regression)
for m = 1:numel(quad1D)
    p_est = polyfit(log(h1), log(err1D_all(m,:)), 1);
    fprintf('%s: p ≈ %.3f\n', quadNames1D{m}, p_est(1));   % p ≈ ordning
end

% Ordning för 2D-trapets
p_est_2D = polyfit(log(h2), log(err2D), 1);
fprintf('Trapets 2D: p ≈ %.3f\n', p_est_2D(1));
fprintf('--------------------------------------\n\n');

%% 4) PLOT AV FELEN
figure; hold on;
for m = 1:numel(quad1D)
    loglog(h1, err1D_all(m,:), markers1D{m}, 'DisplayName', quadNames1D{m});
end
loglog(h2, err2D, 'd-', 'DisplayName','Trapets 2D');

grid on;
xlabel('h');
ylabel('|E_h|');
legend('Location','best');
title('Konvergensordning');
end



%% GENERELLA HJÄLPFUNKTIONER 

function V = volume1D(quadFun, R, alpha, n)
% 2*pi * ∫_0^R [g(R)-g(r)] r dr med vald kvadratur
I = quadFun(@(r) integrand(r,alpha,R), 0, R, n);
V = 2*pi*I;
end

function V = volume2D(alpha,R,L,n)
% Volym via 2D-trapets på kvadrat [-L/2,L/2]^2
h = L/n;
x = linspace(-L/2, L/2, n+1);
y = x;
[X,Y] = meshgrid(x,y);
r = sqrt(X.^2 + Y.^2);

G = g_glass(R,alpha) - g_glass(r,alpha);
G(r>R) = 0;                     % utanför glaset = 0
V = trap2d(G, h, h);
end

function val = integrand(r, alpha, R)
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
if mod(n,2) ~= 0, error('Simpson kräver jämnt n.'); end
x = linspace(a,b,n+1);
h = (b-a)/n;
y = f(x);
I = h/3 * (y(1) + y(end) + 4*sum(y(2:2:end)) + 2*sum(y(3:2:end-1)));
end

function I = trap2d(F,hx,hy)
Irow = hx * (0.5*F(:,1) + sum(F(:,2:end-1),2) + 0.5*F(:,end));
I    = hy * (0.5*Irow(1) + sum(Irow(2:end-1)) + 0.5*Irow(end));
end
