function u4_cNY_flux_Tmed
% Allmän flux + Tmed-beräkning med steglängdshalvering,
% återanvänder solve_bvp och postprocess (BVP-solver-stil).

%% --- PARAMETRAR (ÄNDRA HÄR FÖR ANDRA PROBLEM) ---
a        = 1;      % inre radie
b        = 2;      % yttre radie
T_vatska = 400;    % T(a)
T_luft   = 22;     % T(b)

% ODE: r T'' + T' = 0  <=>  A(r)T'' + B(r)T' + C(r)T = D(r)
A_ode = @(r) r;
B_ode = @(r) 1;
C_ode = @(r) 0;
D_ode = @(r) 0;

tol   = 0.5e-4;    % tolerans för Fin, Fout, Tmed
n0    = 19;        % start-antal inre punkter (som i din gamla kod)

%% --- INITIERA LISTOR FÖR KONVERGENS ---
h_list    = [];
Fin_list  = [];
Fout_list = [];
Tmed_list = [];

% första diskretisering
h = (b-a)/(n0+1);
n = n0;

[r,T] = solve_bvp(h, a, b, T_vatska, T_luft, A_ode, B_ode, C_ode, D_ode);
[Fin,Fout,Tmed] = postprocess(r,T,a,b);

vals_prev = [Fin,Fout,Tmed];

h_list(end+1)    = h;
Fin_list(end+1)  = Fin;
Fout_list(end+1) = Fout;
Tmed_list(end+1) = Tmed;

fprintf('n=%4d, h=%.6f, Fin=%.8f, Fout=%.8f, Tmed=%.8f\n',...
        n,h,Fin,Fout,Tmed);

%% --- STEGLÄNGDSHALVERING (n -> 2n) ---
while true
    n = 2*n;
    h = (b-a)/(n+1);

    [r,T] = solve_bvp(h, a, b, T_vatska, T_luft, A_ode, B_ode, C_ode, D_ode);
    [Fin,Fout,Tmed] = postprocess(r,T,a,b);

    vals_new = [Fin,Fout,Tmed];
    diff = abs(vals_new - vals_prev);
    maxdiff = max(diff);

    h_list(end+1)    = h;
    Fin_list(end+1)  = Fin;
    Fout_list(end+1) = Fout;
    Tmed_list(end+1) = Tmed;

    fprintf('n=%4d, h=%.6f, Fin=%.8f (Δ=%.3e), Fout=%.8f (Δ=%.3e), Tmed=%.8f (Δ=%.3e)\n',...
        n,h,Fin,diff(1),Fout,diff(2),Tmed,diff(3));

    if maxdiff < tol
        break;
    end
    vals_prev = vals_new;
end

%% --- SLUTUTSKRIFT ---
fprintf('\nSlutligt:\nFin=%.8f, Fout=%.8f, skillnad=%.3e\n',Fin,Fout,Fin-Fout);
Tmean_simple = 0.5*(T_vatska + T_luft);
fprintf('Tmed=%.8f, (Tv+Tl)/2=%.8f\n',Tmed,Tmean_simple);

%% --- NOGGRANNHETSORDNING p FÖR Fin, Fout, Tmed ---
if numel(h_list) >= 3
    Fin_ref  = Fin_list(end);
    Fout_ref = Fout_list(end);
    Tmed_ref = Tmed_list(end);

    errFin  = abs(Fin_list  - Fin_ref);
    errFout = abs(Fout_list - Fout_ref);
    errTmed = abs(Tmed_list - Tmed_ref);

    idx = 1:numel(h_list)-1;   % hoppa över sista (err=0)

    pFin  = polyfit(log(h_list(idx)), log(errFin(idx)),  1);
    pFout = polyfit(log(h_list(idx)), log(errFout(idx)), 1);
    pTmed = polyfit(log(h_list(idx)), log(errTmed(idx)), 1);

    fprintf('\nBeräknad noggrannhetsordning (≈ p):\n');
    fprintf('  Fin : p ≈ %.3f\n', pFin(1));
    fprintf('  Fout: p ≈ %.3f\n', pFout(1));
    fprintf('  Tmed: p ≈ %.3f\n', pTmed(1));
end

end

%% === ÅTERANVÄNDA LOKALA FUNKTIONER (samma som i BVP-solvern) ===

function [r, T] = solve_bvp(h, a, b, Ta, Tb, A_ode, B_ode, C_ode, D_ode)
% Löser: A(r)T'' + B(r)T' + C(r)T = D(r)
    n = round((b - a) / h) - 1;
    h = (b - a) / (n + 1); 

    r = a + (0:n+1)' * h;
    r_inner = r(2:end-1);

    Amat = sparse(n,n);
    f = zeros(n,1);

    for k = 1:n
        rk = r_inner(k);

        c_prev = A_ode(rk)/h^2 - B_ode(rk)/(2*h);
        c_curr = -2*A_ode(rk)/h^2 + C_ode(rk);
        c_next = A_ode(rk)/h^2 + B_ode(rk)/(2*h);
        f(k)   = D_ode(rk);

        Amat(k,k) = c_curr;

        if k > 1
            Amat(k,k-1) = c_prev;
        else
            f(k) = f(k) - c_prev*Ta;
        end

        if k < n
            Amat(k,k+1) = c_next;
        else
            f(k) = f(k) - c_next*Tb;
        end
    end

    T_inner = Amat \ f;
    T = [Ta; T_inner; Tb];
end

function [F_in, F_out, T_medel] = postprocess(r, T, a, b)
    h = r(2) - r(1);

    % derivator vid r=a och r=b (2:a ordningens 3-punktsformler)
    dT_a = (-3*T(1)   + 4*T(2)     - T(3))      / (2*h);
    dT_b = ( 3*T(end) - 4*T(end-1) + T(end-2))  / (2*h);

    F_in  = -2*pi*a*dT_a;
    F_out = -2*pi*b*dT_b;

    % medeltemperatur (vanlig T_med, ej kvadratmedel):
    integrand = T .* r;
    integral_val = trapz(r, integrand);
    T_medel = 2/(b^2 - a^2) * integral_val;
end
