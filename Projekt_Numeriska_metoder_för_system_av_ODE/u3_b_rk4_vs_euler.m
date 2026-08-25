function u3_b_rk4_vs_euler
% Uppgift 3(b): Framåt Euler och RK4, bestäm x(20), y(20)
% med fel < 1e-2 och rapportera steglängder. RK4-approx med
% mycket fin h används som "referenslösning".

eta  = 0.05;
beta = 0.10;
u0   = [-1; 0; -3; 0];
t0   = 0;
T    = 20;

f = @f_springs;

%% 1. Referenslösning med fin RK4
h_ref = 0.0005;   % väldigt fin steglängd
[~,Uref] = rk4_system(f,t0,T,u0,h_ref,eta,beta);
xRef = Uref(1,end); 
yRef = Uref(3,end);

fprintf("Referens (RK4, h=%.4g): x(20)=%.6f, y(20)=%.6f\n\n", ...
        h_ref, xRef, yRef);

tol = 1e-2;
maxIter = 30;

%% 2. Euler: öka N (fördubbla) tills fel < 1e-2
fprintf("=== Euler: N-fördubbling (h = (T-t0)/N) ===\n");

N_E      = 1000;   % start med 1000 steg
best_h_E = NaN;
best_N_E = NaN;

% ---- NYTT: lagra h och fel för Euler ----
h_list_E   = [];
err_list_E = [];

for it = 1:maxIter
    hE = (T - t0)/N_E;   % steglängd h baserad på N
    [~,U] = euler_system(f,t0,T,u0,hE,eta,beta);
    x = U(1,end); 
    y = U(3,end);
    err = max(abs([x - xRef, y - yRef]));

    fprintf("Euler: N = %6d   h = %.6f   x(20)=%.6f  y(20)=%.6f  max-fel=%.3e\n", ...
            N_E, hE, x, y, err);

    % ---- NYTT: spara h och fel ----
    h_list_E(end+1)   = hE;
    err_list_E(end+1) = err;

    if err < tol
        best_h_E = hE;
        best_N_E = N_E;
        break;
    end

    N_E = 2*N_E;  % fördubbla antal steg
end

if ~isnan(best_h_E)
    fprintf("Euler: minsta N (med fördubbling) som ger fel < 1e-2: N = %d, h ≈ %.6f\n\n", ...
            best_N_E, best_h_E);
else
    fprintf("Euler: nådde inte fel < 1e-2 inom %d fördubblingar av N\n\n", maxIter);
end

% ---- NYTT: beräkna noggrannhetsordningen för Euler ----
if numel(h_list_E) >= 2
    p_E = polyfit(log(h_list_E), log(err_list_E), 1);
    fprintf("Noggrannhetsordning (Euler): p ≈ %.3f\n\n", p_E(1));
end


%% 3. RK4: öka N (fördubbla) tills fel < 1e-2
fprintf("=== RK4: N-fördubbling (h = (T-t0)/N) ===\n");

N_R      = 50;   % start med 1000 steg
best_h_R = NaN;
best_N_R = NaN;

% ---- NYTT: lagra h och fel för RK4 ----
h_list_R   = [];
err_list_R = [];

for it = 1:maxIter
    hR = (T - t0)/N_R;   % steglängd h baserad på N
    [~,U] = rk4_system(f,t0,T,u0,hR,eta,beta);
    x = U(1,end); 
    y = U(3,end);
    err = max(abs([x - xRef, y - yRef]));

    fprintf("RK4  : N = %6d   h = %.6f   x(20)=%.6f  y(20)=%.6f  max-fel=%.3e\n", ...
            N_R, hR, x, y, err);

    % ---- NYTT: spara h och fel ----
    h_list_R(end+1)   = hR;
    err_list_R(end+1) = err;

    if err < tol
        best_h_R = hR;
        best_N_R = N_R;
        break;
    end

    N_R = 2*N_R;  % fördubbla antal steg
end

if ~isnan(best_h_R)
    fprintf("RK4 : minsta N (med fördubbling) som ger fel < 1e-2: N = %d, h ≈ %.6f\n\n", ...
            best_N_R, best_h_R);
else
    fprintf("RK4 : nådde inte fel < 1e-2 inom %d fördubblingar av N\n\n", maxIter);
end

% ---- NYTT: beräkna noggrannhetsordningen för RK4 ----
if numel(h_list_R) >= 2
    p_R = polyfit(log(h_list_R), log(err_list_R), 1);
    fprintf("Noggrannhetsordning (RK4): p ≈ %.3f\n\n", p_R(1));
end


%% 4. Plotta en RK4-lösning (t.ex. med det h som klarade kravet)
if ~isnan(best_h_R)
    hPlot = best_h_R;
else
    hPlot = 0.01;  % fallback
end

[tPlot,UPlot] = rk4_system(f,t0,T,u0,hPlot,eta,beta);
xPlot = UPlot(1,:); 
yPlot = UPlot(3,:);

figure;
plot(tPlot,xPlot,'b',tPlot,yPlot,'r'); grid on;
xlabel('t'); ylabel('position');
legend('x(t)','y(t)');
title(sprintf('RK4-lösning, h=%.4f',hPlot));

makeVideo = false;   % true för video
if makeVideo
    animatesprings(tPlot, xPlot, yPlot, 3, 'springs.mp4');
end
end

%% ===== HJÄLPFUNKTIONER =====
function dudt = f_springs(t,u,eta,beta)
x = u(1); v = u(2); y = u(3); w = u(4);
dxdt = v;
dvdt = -eta*v - 2*x - beta*x^3 + y - 1;
dydt = w;
dwdt = -eta*w + 0.5*(x - y) - 1;
dudt = [dxdt; dvdt; dydt; dwdt];
end

function [t,U] = euler_system(f,t0,T,u0,h,eta,beta)
N = round((T - t0)/h);
h = (T - t0)/N;           % säkerställ att T/h blir exakt
t = linspace(t0,T,N+1);
U = zeros(length(u0),N+1);
U(:,1) = u0;
for k=1:N
    U(:,k+1) = U(:,k) + h * f(t(k),U(:,k),eta,beta);
end
end

function [t,U] = rk4_system(f,t0,T,u0,h,eta,beta)
N = round((T - t0)/h);
h = (T - t0)/N;           % säkerställ att T/h blir exakt
t = linspace(t0,T,N+1);
U = zeros(length(u0),N+1);
U(:,1) = u0;
for k=1:N
    k1 = f(t(k),    U(:,k),          eta,beta);
    k2 = f(t(k)+h/2,U(:,k)+h/2*k1,   eta,beta);
    k3 = f(t(k)+h/2,U(:,k)+h/2*k2,   eta,beta);
    k4 = f(t(k)+h,  U(:,k)+h*k3,     eta,beta);
    U(:,k+1) = U(:,k) + h*(k1+2*k2+2*k3+k4)/6;
end
end



