function u1_dNY_find_alpha
% Beräknar alpha både för rotationssymmetriska glaset och glaset med kvadratisk öppning

%% Probleminställningar
R       = 3;              % radie
L       = 3*sqrt(2);      % kvadratens sida
Vtarget = 100;            % önskad volym

nInt1D  = 400;           % noggrannhet rotation
nInt2D  = 200;           % noggrannhet 2D kvadrat

% Startgissningar för sekantmetod (samma för båda)
alpha0  = 0.6;
alpha1  = 0.9;

tolF    = 1e-6;
tolX    = 1e-6;
maxIt   = 50;


%% ======= ROTATIONSSYMMETRISKT GLAS =======
fprintf("\n=== Rotationssymmetriskt glas ===\n");

F_rot = @(a) volume_rotation(a, R, nInt1D) - Vtarget;
[alpha_rot, it_rot] = secant(F_rot, alpha0, alpha1, tolF, tolX, maxIt);
V_rot = volume_rotation(alpha_rot, R, nInt1D);

fprintf("Alpha (rotation) ≈ %.10f   (iter = %d)\n", alpha_rot, it_rot);
fprintf("Kontrollvolym ≈ %.10f   |F(alpha)| ≈ %.3e\n\n", V_rot, F_rot(alpha_rot));


%% ======= KVADRATISKT GLAS =======
fprintf("\n=== Kvadratisk öppning ===\n");

F_sq = @(a) volume_square(a, R, L, nInt2D) - Vtarget;
[alpha_sq, it_sq] = secant(F_sq, alpha0, alpha1, tolF, tolX, maxIt);
V_sq = volume_square(alpha_sq, R, L, nInt2D);

fprintf("Alpha (kvadrat) ≈ %.10f   (iter = %d)\n", alpha_sq, it_sq);
fprintf("Kontrollvolym ≈ %.10f   |F(alpha)| ≈ %.3e\n\n", V_sq, F_sq(alpha_sq));

end


%% ======== FUNKTIONER =========

function V = volume_rotation(alpha, R, n)
% 2π∫0^R (g(R,α) - g(r,α))·r dr
a = 0; b = R;
f = @(r) (g_glass(R,alpha) - g_glass(r,alpha)) .* r;
I = simpson1d(f, a, b, n);
V = 2*pi*I;
end

function V = volume_square(alpha, R, L, n)
% Volym via 2D-trapets i [-L/2,L/2]^2
h = L/n;
x = linspace(-L/2, L/2, n+1);
y = x;
[X,Y] = meshgrid(x,y);
r = sqrt(X.^2 + Y.^2);

G = g_glass(R,alpha) - g_glass(r,alpha);
G(r > R) = 0;

V = trap2d(G, h, h);
end

function g = g_glass(r, alpha)
R = 3;
g = 6 .* r.^3 .* exp(-r) .* (3 + 2*sin(alpha*R)) ./ (3 + 2*sin(alpha.*r));
end


%% ===== Numerisk integration =====

function I = simpson1d(f, a, b, n)
if mod(n,2) ~= 0, n = n+1; end % fixar ojämnt n
x = linspace(a,b,n+1);
h = (b-a)/n;
y = f(x);
I = h/3 * (y(1) + y(end) + 4*sum(y(2:2:end)) + 2*sum(y(3:2:end-1)));
end

function I = trap2d(F, hx, hy)
Irow = hx * (0.5*F(:,1) + sum(F(:,2:end-1),2) + 0.5*F(:,end));
I    = hy * (0.5*Irow(1) + sum(Irow(2:end-1)) + 0.5*Irow(end));
end


%% ====== Sekantmetod ======
function [x2, it] = secant(F, x0, x1, tolF, tolX, maxIt)
f0 = F(x0);
f1 = F(x1);

for it = 1:maxIt
    if abs(f1 - f0) < eps, x2 = x1; return; end
    x2 = x1 - f1*(x1 - x0)/(f1 - f0);
    f2 = F(x2);
    if abs(f2) < tolF || abs(x2 - x1) < tolX, return; end
    x0 = x1; f0 = f1;
    x1 = x2; f1 = f2;
end
warning('Sekant: max iterationer nådda.');
end