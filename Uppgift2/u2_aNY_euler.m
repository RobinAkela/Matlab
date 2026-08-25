function u2_aNY_euler
%% PROBLEM-PARAMETRAR (ÄNDRA HÄR FÖR ANDRA ODE:er)
f      = @(t,y) sin(pi*t) - 2*y + (4/3)*y.^2;  % ODE
y0     = 1.2;                                    % begynnelsevärde
t0     = 0;                                      % starttid
T      = 8;                                      % sluttid (vi vill ha y(T))
tol    = 1e-4;                                   % felkrav på y(T)
h0     = 0.32;                                    % startsteg (ger 4 steg här)
maxHalv = 20;                                    % max antal halveringar
doPlot = true;                                  % sätt false om du inte vill plotta

%% START: n från h0 
n0 = round((T - t0)/h0);       % antal steg för startsteget
n  = n0;

[t,y] = euler_forward(f,t0,T,y0,n);
y_old = y(end);

fprintf("Start: n = %d, h = %.5f, y(T) ≈ %.10f\n", n, (T-t0)/n, y_old);

% lagra för plott
t_list{1} = t;
y_list{1} = y;
n_list(1) = n;

%% STEGLÄNGDSHALVERING 
halvningar = 0;
errEst = inf;

while errEst > tol && halvningar < maxHalv
    halvningar = halvningar + 1;
    n = 2*n;                        % halverar h → dubblar n
    h = (T - t0)/n;

    [t,y] = euler_forward(f,t0,T,y0,n);
    y_new = y(end);

    errEst = abs(y_new - y_old);

    fprintf("Halvning %2d: n = %5d, h = %.6f, y(T) ≈ %.10f, |Δ| = %.3e\n", ...
        halvningar, n, h, y_new, errEst);

    % spara för plott
    t_list{end+1} = t;
    y_list{end+1} = y;
    n_list(end+1) = n;

    y_old = y_new;
end

fprintf("\nAntal halveringar som krävs: %d\n", halvningar);
fprintf("Slutapproximation av y(%g) ≈ %.10f\n\n", T, y_new);

%%   PLOTT AV LÖSNINGAR FÖR ALLA n  
if doPlot
    figure; hold on;
    for k = 1:numel(n_list)
        plot(t_list{k}, y_list{k}, 'DisplayName', sprintf('n=%d', n_list(k)));
    end
    xlabel('t'); ylabel('y(t)');
    title('Framåt Euler – lösning för olika n (steglängdshalvering)');
    legend('Location','best');
    grid on;
end
end

%%   FRAMÅT EULER  
function [t,y] = euler_forward(f,t0,T,y0,n)
h = (T - t0)/n;
t = linspace(t0,T,n+1);
y = zeros(size(t));
y(1) = y0;
for k = 1:n
    y(k+1) = y(k) + h*f(t(k),y(k));
end
end
