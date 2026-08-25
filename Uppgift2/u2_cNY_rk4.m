function u2_cNY_rk4

%%   Problemparametrar  
f       = @(t,y) sin(pi*t) - 2*y + (4/3)*y.^2;
T       = 8;
t0      = 0;
y0      = 1.2;

tol     = 1e-6;     % felkrav
h       = 0.32;     % startsteg
maxSteps = 20;      % säkerhetsgräns

%%   Första körningen  
n = round((T - t0)/h);
[~,y] = rk4_solver(f,t0,T,y0,n);

y_vals = y(end);
h_vals = h;
n_vals = n;

fprintf("  RK4 med systematisk halvering  \n");
fprintf("n      h           y(T)                 est.err\n");
fprintf("%5d  %10.6f  %.12f      ---\n", n, h, y_vals(end));

%%   Halveringsloop  
errEst = inf;

while errEst > tol && length(h_vals) < maxSteps
    % halvera h → dubblera n
    h = h/2;
    n = round((T - t0)/h);

    [~,y_new] = rk4_solver(f,t0,T,y0,n);
    
    errEst = abs(y_new(end) - y_vals(end));

    fprintf('%5d  %10.6f  %.12f   %.3e\n', n, h, y_new(end), errEst);

    % spara historik
    y_vals(end+1) = y_new(end);
    h_vals(end+1) = h;
    n_vals(end+1) = n;
end

%%   Resultat  
fprintf("\nStörsta h som uppfyller tol < %.1e: h = %.6f (n = %d)\n", ...
         tol, h_vals(end), n_vals(end));

%%   Konvergensordning (q & p)  
if length(y_vals) >= 3
    fprintf("\nKvoter q och ordning p ≈ log2|q|:\n");
    for k = 3:length(y_vals)
        Y4h = y_vals(k-2);
        Y2h = y_vals(k-1);
        Yh  = y_vals(k);

        q = (Y4h - Y2h) / (Y2h - Yh);
        p = log2(abs(q));
        fprintf('[%6d,%6d,%6d] → q = %.5f,  p ≈ %.3f\n', ...
                n_vals(k-2), n_vals(k-1), n_vals(k), q, p);
    end
end

end

%%   Hjälpfunktioner  
function [t,y] = rk4_solver(f,t0,T,y0,n)
h = (T-t0)/n;
t = linspace(t0,T,n+1);
y = zeros(size(t)); 
y(1)=y0;

for k=1:n
    tk = t(k); 
    yk = y(k);
    
    k1 = f(tk, yk);
    k2 = f(tk + h/2, yk + h*k1/2);
    k3 = f(tk + h/2, yk + h*k2/2);
    k4 = f(tk + h,   yk + h*k3);

    y(k+1) = yk + h*(k1 + 2*k2 + 2*k3 + k4)/6;
end
end
