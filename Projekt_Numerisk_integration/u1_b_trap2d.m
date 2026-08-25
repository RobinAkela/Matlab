function u1_b_trap2d
% 2D-trapetsregel för glasvolym med kvadratisk öppning + steglängdshalvering

%% PARAMETRAR
alpha  = 0.9;           % formfaktor
R      = 3;             % glasradie

L      = 3*sqrt(2);     % sidlängd på kvadraten
xmin   = -L/2;          % x-intervall
xmax   =  L/2;
ymin   = -L/2;          % y-intervall
ymax   =  L/2;

nStart = 30;            % start-n (antal delintervall i x- och y-led)
tol    = 0.5e-3;        % mål: 3 korrekta decimaler i volymen

%% STARTVÄRDE
volFun = @(n) volume2D(alpha,R,xmin,xmax,ymin,ymax,n);

n     = nStart;
V_old = volFun(n);
hx    = (xmax-xmin)/n;

fprintf("n = %4d, h = %.6f, V ≈ %.10f\n", n, hx, V_old);

%% STEGLÄNGDSHALVERING 
while true
    n     = 2*n;              % halvera h i både x- och y-led
    V_new = volFun(n);
    hx    = (xmax-xmin)/n;
    diff  = abs(V_new - V_old);

    fprintf("n = %4d, h = %.6f, V ≈ %.10f, |Δ| ≈ %.3e\n", ...
            n, hx, V_new, diff);

    if diff < tol
        break;
    end
    V_old = V_new;
end

fprintf("\nSlutlig volym (2D trapets, 3 dec): V ≈ %.10f (n = %d)\n", V_new, n);
end

%% GENERELLA HJÄLPFUNKTIONER
function V = volume2D(alpha,R,xmin,xmax,ymin,ymax,n)
% Beräknar volymen för givna alpha, R och kvadrat [xmin,xmax]×[ymin,ymax]

hx = (xmax-xmin)/n;
hy = (ymax-ymin)/n;

x = linspace(xmin, xmax, n+1);
y = linspace(ymin, ymax, n+1);
[X,Y] = meshgrid(x,y);

r = sqrt(X.^2 + Y.^2);                 % avstånd till centrum

G = g_glass(R,alpha) - g_glass(r,alpha);
G(r > R) = 0;                          % sätt allt utanför glaset till 0

V = trap2d(G, hx, hy);
end

function g = g_glass(r, alpha)
R = 3;
g = 6 .* r.^3 .* exp(-r) .* (3 + 2*sin(alpha*R)) ./ (3 + 2*sin(alpha*r));
end

function I = trap2d(F, hx, hy)
% 2D-trapets: först i x-led radvis, sedan i y-led
Irow = hx * (0.5*F(:,1) + sum(F(:,2:end-1),2) + 0.5*F(:,end));
I    = hy * (0.5*Irow(1) + sum(Irow(2:end-1)) + 0.5*Irow(end));
end
