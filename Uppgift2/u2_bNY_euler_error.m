function u2_bNY_euler_error

%%   Problemparametrar  
f      = @(t,y) sin(pi*t) - 2*y + (4/3)*y.^2;
T      = 8;
t0     = 0;
y0     = 1.2;

tol    = 1e-3;     % felkrav
h      = 0.32;      % startsteg (kan ändras)
maxSteps = 20;     % säkerhetsgräns

%%   Första lösningen  
n = round((T - t0)/h);
[~,y] = euler_forward(f,t0,T,y0,n);

y_vals = y(end);
h_vals = h;
n_vals = n;

fprintf("n      h        y(T)          est.err\n");
fprintf("%5d  %7.5f  %.10f      ---\n", n, h, y(end));

%%   Halveringsloop  
errEst = inf;

while errEst > tol && length(h_vals) < maxSteps
    
    % halvera h -> dubblera n
    h = h/2;
    n = round((T - t0)/h);
    [~,y_new] = euler_forward(f,t0,T,y0,n);

    errEst = abs(y_new(end) - y_vals(end));

    fprintf('%5d  %7.5f  %.10f  %.3e\n', n, h, y_new(end), errEst);

    % spara historik
    y_vals(end+1) = y_new(end);
    h_vals(end+1) = h;
    n_vals(end+1) = n;
end

%%   Resultat  
fprintf("\nStörsta h med fel < %.1e: h = %.5f (n = %d)\n", tol, h_vals(end), n_vals(end));

%%   Ordning (kvoter)  
if length(y_vals) >= 3
    fprintf("\nKvoter q och p ≈ log2|q|:\n");
    for k = 3:length(y_vals)
        Y4h = y_vals(k-2);
        Y2h = y_vals(k-1);
        Yh  = y_vals(k);
        q = (Y4h-Y2h)/(Y2h-Yh);
        p = log2(abs(q));
        fprintf('[%4d,%4d,%4d]  q=%.4f, p=%.3f\n', n_vals(k-2),n_vals(k-1),n_vals(k), q, p);
    end
end

end


%%   Hjälpfunktioner  
function dydt = f_ode(t,y)
dydt = sin(pi*t) - 2*y + (4/3)*y.^2;
end

function [t,y] = euler_forward(f,t0,T,y0,n)
h = (T-t0)/n;
t = linspace(t0,T,n+1);
y = zeros(size(t)); 
y(1)=y0;
for k=1:n
    y(k+1) = y(k) + h*f(t(k),y(k));
end
end
