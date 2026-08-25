% BVP Solver 
% Löser: r * T'' + T' = 0, T(a) = T_vätska, T(b) = T_luft

clear; clc; close all;

%% --- PARAMETRAR (ÄNDRA HÄR) ---
a = 1;                  % Inre radie
b = 2;                  % Yttre radie
T_vatska = 400;         % Randvillkor vid r=a
T_luft = 22;            % Randvillkor vid r=b

h_start = 0.2;         % Start-steglängd 
% n_fix = 4;           % om ger n direkt (då ignoreras h_start)


A_ode = @(r) r;         % Koefficient framför T''
B_ode = @(r) 1;         % Koefficient framför T'
C_ode = @(r) 0;         % Koefficient framför T
D_ode = @(r) 0;         % Högerled

% Extra konvergensloop (stäng av/på här)
KONVERGENS_TILL_4_SIFFROR = true;  % Sätt till true för automatisk konvergens till 4 korrekta siffror

%% --- BERÄKNINGAR ---
% Om n_fix är definierad, beräkna h_start från n istället
if exist('n_fix', 'var')
    h_start = (b - a) / (n_fix + 1);
end

% Skapa tre steglängder: h, h/2, h/4
h_vec = [h_start, h_start/2, h_start/4]; % vill göra samma beräkning för tre olika steglängder
% syfte kunna jmf resultat för olika h sedan uppskatta konvergensordningen

fprintf('%-10s %-8s %-15s %-15s %-15s %-15s\n', 'h', 'n', 'F_in', 'F_out', '|F_in-F_out|', 'T_medel');
fprintf('%s\n', repmat('-', 1, 80));

F_in_values = [];

if KONVERGENS_TILL_4_SIFFROR
    % === AUTOMATISK KONVERGENS TILL 4 KORREKTA SIFFROR ===
    % Fortsätter halvera h tills felet < 0.5e-4 (minst 4 korrekta siffror)
    % Feluppskattning: |V_h - V| ≤ |V_{2h} - V_h|
    
    tol = 0.5e-4;  % Tolerans för 4 korrekta siffror
    h_current = h_start;
    F_in_prev = [];
    F_out_prev = [];
    T_medel_prev = [];
    
    iter = 0;
    max_iter = 15;  % Säkerhetsgräns
    
    while iter < max_iter
        iter = iter + 1;
        n = round((b - a) / h_current) - 1;
        h_current = (b - a) / (n + 1);  % Justera h för exakt n
        
        % Lös BVP
        [r, T] = solve_bvp(h_current, a, b, T_vatska, T_luft, A_ode, B_ode, C_ode, D_ode);
        
        % Beräkna storheter
        [F_in, F_out, T_medel] = postprocess(r, T, a, b);
        F_in_values(end+1) = F_in;
        
        % Skriv ut resultat
        fprintf('%-10.6f %-8d %-15.6f %-15.6f %-15.6e %-15.6f\n', ...
            h_current, n, F_in, F_out, abs(F_in - F_out), T_medel);
        
        % Kontrollera konvergens (feluppskattning: |V_h - V| ≤ |V_{2h} - V_h|)
        if ~isempty(F_in_prev)
            err_F_in = abs(F_in - F_in_prev);
            err_F_out = abs(F_out - F_out_prev);
            err_T_medel = abs(T_medel - T_medel_prev);
            
            % Om alla fel är mindre än toleransen, är vi klara
            if err_F_in < tol && err_F_out < tol && err_T_medel < tol
                fprintf('\nKonvergens uppnådd! Fel < %.1e för alla storheter.\n', tol);
                fprintf('Fel: F_in=%.2e, F_out=%.2e, T_medel=%.2e\n', ...
                    err_F_in, err_F_out, err_T_medel);
                break;
            end
        end
        
        % Spara värden för nästa iteration
        F_in_prev = F_in;
        F_out_prev = F_out;
        T_medel_prev = T_medel;
        
        % Halvera h för nästa iteration
        h_current = h_current / 2;
    end
    
    if iter >= max_iter
        fprintf('\nVarning: Max antal iterationer nått. Konvergens kanske inte uppnådd.\n');
    end
    
    % Spara h och n för plottning
    h_plot = h_current;
    n_plot = n;
    
else
    % === ORIGINAL BERÄKNING MED TRE STEGLÄNGDER ===
        for i = 1:length(h_vec) %aktuellt h värde från vektorn (3)
        h = h_vec(i);
        n = round((b - a) / h) - 1;
        h = (b - a) / (n + 1);  % Justera h för exakt n, KONSEKVENS: h i utskriften är inte exakt start värdet korresponderar till n så nära som möjligt
        
        % Lös BVP anropa lösningsfkn: bygger matris, högerled f
        [r, T] = solve_bvp(h, a, b, T_vatska, T_luft, A_ode, B_ode, C_ode, D_ode); 
        
        % Beräkna storheter (använder resultatet från förra raden) 
        [F_in, F_out, T_medel] = postprocess(r, T, a, b); %fkn postprocess approximerar T'(a) och T'(b), beräknarnar Fin/Fout och Tmedel med trapetsregeln
        F_in_values(i) = F_in; % sparar F_in för varje h (behövs när vi ska beräkna empiriska ordningen p)
        
        % Skriv ut resultat
        fprintf('%-10.6f %-8d %-15.6f %-15.6f %-15.6e %-15.6f\n', ...
            h, n, F_in, F_out, abs(F_in - F_out), T_medel);
    end
    
    % Spara h och n för plottning (sista iterationen)
    h_plot = h;
    n_plot = n;
end



if length(F_in_values) >= 3 % för att kunna uppskatta konv.ord. p krävs minst 3 storlekar  ((KAN DET VARA MER HUR LÖSER MAN DET DÅ?
    diff1 = abs(F_in_values(1) - F_in_values(2)); % ungefärligt fel för grövre g
    diff2 = abs(F_in_values(2) - F_in_values(3)); % ungefärligt fel för det halverade h
    if diff2 > 0
        p = log2(diff1 / diff2);
        fprintf('\nEmpirisk ordning p = %.2f\n', p);
    end
end

%%Jämför med (T_vätska + T_luft)/2 för finaste h
%OBS i en cylindrisk geometri väger den kalla sidan mer eftersom arean ökar med r
T_medel_simple = (T_vatska + T_luft)/2;
fprintf('\n(T_vätska + T_luft)/2 = %.6f\n', T_medel_simple);
fprintf('T_medel (finaste h)   = %.6f\n', T_medel);
if T_medel > T_medel_simple
    fprintf('T_medel är större än (T_vätska + T_luft)/2.\n'); % kan inte stämma för denna geometri
else
    fprintf('T_medel är mindre än (T_vätska + T_luft)/2.\n');
end

% Enkel kontroll av randvillkor och att T avtar (finaste h)
fprintf('\nKontroll för finaste h:\n');
fprintf('T(a) = %.6f (förväntat %.6f)\n', T(1), T_vatska);
fprintf('T(b) = %.6f (förväntat %.6f)\n', T(end), T_luft);
if all(diff(T) < 0) % diff T beräknar skillnaden melan varje par av efterföljande element i vektorn T
    fprintf('Temperaturen avtar monotoniskt med r.\n'); %om diffen mellan alla är negativ -> avtar monotoniskt
else
    fprintf('Temperaturen är inte strikt avtagande överallt.\n');
end 

%% --- PLOTT ---
% r och T kommer från sista beräkningen (finaste h)
figure;
plot(r, T, '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
hold on;
plot(a, T_vatska, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot(b, T_luft, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('r');
ylabel('T(r)');
title(sprintf('Temperaturfördelning (h=%.4f, n=%d)', h_plot, n_plot));
grid on;
legend('T(r)', 'Randvillkor', 'Location', 'best');

%% Utskrift av T vid specifika r-punkter
r_ut = [1.0, 1.25, 1.5, 1.75, 2.0];
fprintf('\nT vid specifika punkter:\n');
for i = 1:length(r_ut)
    [~, idx] = min(abs(r - r_ut(i)));
    fprintf('T(%.2f) = %.6f\n', r(idx), T(idx));
end

%% Utskrift av F (värmeflöde) vid specifika r-punkter
% F(r) = -2*pi*r*T'(r) för cylindergeometri
% ÄNDRA r_ut_F för andra punkter på provet
r_ut_F = [1.0, 1.5, 2.0];  
h_grid = r(2) - r(1);
fprintf('\nF (värmeflöde) vid specifika punkter:\n');
for i = 1:length(r_ut_F)
    [~, idx] = min(abs(r - r_ut_F(i)));
    
    % Beräkna T'(r) med lämplig differensformel
    if idx == 1
        % Vänster rand: framåtdifferens
        dT = (-3*T(idx) + 4*T(idx+1) - T(idx+2)) / (2*h_grid);
    elseif idx == length(r)
        % Höger rand: bakåtdifferens
        dT = (3*T(idx) - 4*T(idx-1) + T(idx-2)) / (2*h_grid);
    else
        % Inre punkt: centraldifferens
        dT = (T(idx+1) - T(idx-1)) / (2*h_grid);
    end
    
    F_r = -2*pi*r(idx)*dT;
    fprintf('F(%.2f) = %.6f\n', r(idx), F_r);
end

%% --- LOKALA FUNKTIONER ---

function [r, T] = solve_bvp(h, a, b, Ta, Tb, A_ode, B_ode, C_ode, D_ode)
% Löser ODE: A(r)*T'' + B(r)*T' + C(r)*T = D(r) med centraldifferenser
% h - steglängd, a,b - ändpunkter, Ta,Tb - randvärden
% A_ode, B_ode, C_ode, D_ode - funktionshandtag för ODE-koefficienterna
   
    n = round((b - a) / h) - 1;
    h = (b - a) / (n + 1); 

    % Bygg r-vektorn
    r = a + (0:n+1)' * h;
    r_inner = r(2:end-1);
    
    % Bygg matris Amat och högerled f
    Amat = sparse(n, n);
    f = zeros(n, 1);
    
    for k = 1:n
        rk = r_inner(k);

        c_prev = A_ode(rk)/h^2 - B_ode(rk)/(2*h);   % Koeff för T_{k-1}
        c_curr = -2*A_ode(rk)/h^2 + C_ode(rk);       % Koeff för T_k
        c_next = A_ode(rk)/h^2 + B_ode(rk)/(2*h);   % Koeff för T_{k+1}
        f(k) = D_ode(rk);                            % Högerled
        
        Amat(k, k) = c_curr;
        
        if k > 1
            Amat(k, k-1) = c_prev;
        else
            f(k) = f(k) - c_prev * Ta;  % Vänster randvillkor
        end
        
        if k < n
            Amat(k, k+1) = c_next;
        else
            f(k) = f(k) - c_next * Tb;  % Höger randvillkor
        end
    end
   
    % Lös system
    T_inner = Amat \ f;
    T = [Ta; T_inner; Tb];
end

function [F_in, F_out, T_medel] = postprocess(r, T, a, b)
    % Beräknar värmeflöde och medeltemperatur
    
    h = r(2) - r(1); 
    
    % 2:a ordningens 3-punktsformler för derivator
    % T′(x0​)≈(−3T0​+4T1​−T2)/2h (frammåt)​
    dT_a = (-3*T(1) + 4*T(2) - T(3)) / (2*h); %vänstra randpunkt frammåtdifferens
   
    % T′(xn+1​)≈(3Tn+1​−4Tn​+Tn−1)/2h​​
    dT_b = (3*T(end) - 4*T(end-1) + T(end-2)) / (2*h); %högra randpunkt bakåtdifferens
    
    F_in = -2 * pi * a * dT_a; % total värme in genom r=a
    F_out = -2 * pi * b * dT_b; % total värme ut genom r = b
    
    % Medeltemperaturen med trapetsregeln (2:a ordningen) Lägre än Tv+Tl /
    % 2 pga större volym i cirkeln mot ytan

    integrand = T .* r; % .* - elementvis multiplikation det man integrerar
    integral_val = trapz(r, integrand); % nummeriska approx av integralen: T(r)rdr
    T_medel = (2 / (b^2 - a^2)) * integral_val; % gör om integralen till ett korrekt medelvärde enligt cylindergeometrins area viktning
end