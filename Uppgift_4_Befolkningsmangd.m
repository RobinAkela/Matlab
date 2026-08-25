load swedpop.mat

% Vektorerna T (år sedan 1920) och Y (Miljoner)
t = T(:);
y = Y(:);
n = numel(t); % Antal datapunkter

% Linjär modell y(t) = c0 + c1*t => A*c ≈ y
A = [ones(n,1), t];

% Minstakvadratslösning
c = A\y;     % QR-Baserad leastsquares med backslash
c0 = c(1);   % Intercept (miljoner)
c1 = c(2);   % Lutning c1 i miljoner/år

% Anpassade värden, överblivande och norm
yhat = A*c;         % Modellens förutsägelse på alla t_i
r = y - yhat;       % residual r_i = y_i - ŷ_i.
R2norm = norm(r,2); % ||R||_2, mäter total felstorlek

fprintf('c0 = %.6f\nc1 = %.6f miljoner/år\n||R||_2 = %.6f\n', c0, c1, R2norm);

% Figur
figure
plot(1920+t, y, 'o');                     % Plotta givna datapunkter
hold on;
plot(1920+t, yhat, '-', 'LineWidth',1.5); % Plotta linjära modellen
grid on;
xlabel('År');
ylabel('Miljoner invånare');
legend('Data', sprintf('y(t)=%.f + %.3f t',c0, c1), 'Location', 'northwest');
title('linjär minstakvadratanpassning');

% Delfråga b)

% Plotta residualen för linjära modellen
figure
plot(1920+t, r, '-o'); grid on
xlabel('År'); ylabel('Residual (M)')
title('Residual för linjär modell')

% Uppskatta period L automatiskt via FFT (årsdata: Δt≈1)
dt = mean(diff(t));              % Samplingsintervall i år
R  = r - mean(r);                % Tar bort medelvärde (DC)
N  = numel(R);                   % Antal prov
F  = fft(R);                     % Diskret Fouriertransform
freq = (0:N-1)/(N*dt);           % Frekvenser i cykler/år
amp  = abs(F)/N;                 % Amplitudspektrum

halfIdx = 2:floor(N/2);          % hoppa över DC-komponenten
[~, j]  = max(amp(halfIdx));     % Välj dominerande frekvens
f0 = freq(halfIdx(j));           % Huvudfrekvens f0 i cykler/år
L  = 1/f0;                       % period L i år
k  = 2*pi/L;                     % Omvandling till vinkelhastighet k (rad/år)

% Förfinad modell: y = d0 + d1*t + d2*sin(k t) + d3*cos(k t)
B = [ones(n,1),  t,  sin(k*t),  cos(k*t)]; % Ny designmatris med sin/cos på k*t.
d = B\y;                  % Least squares för koefficienterna d0...d3.    
d0 = d(1); 
d1 = d(2); 
d2 = d(3); 
d3 = d(4);

yhat2  = B*d;         % Nya modellens förutsägelse
r2     = y - yhat2;   % Nya residualer
R2norm2 = norm(r2,2); % Jämförelsemål

% Utskfrift
fprintf('\nFörfinad modell:\n')
fprintf('d0 = %.6f\nd1 = %.6f miljoner/år\n', d0, d1);
fprintf('d2 = %.6f\nd3 = %.6f\n', d2, d3);
fprintf('L ≈ %.2f år, k = %.6f 1/år\n', L, k);
fprintf('||R_lin||_2 = %.6f,  ||R_ref||_2 = %.6f\n', R2norm, R2norm2);

% Plotta data + förfinad kurva
figure
plot(1920+t, y, 'o');
hold on
grid on
plot(1920+t, yhat2, '-', 'LineWidth', 1.5);
xlabel('År'); ylabel('Miljoner invånare')
legend('Data', 'Förfinad modell', 'Location','northwest')
title('d0 + d1 t + d2 sin(k t) + d3 cos(k t)')

% Jämför residualer (Visar förbättringen)
figure
plot(1920+t, r, '-o'); 
plot(1920+t, r2, '-'); 
hold on
grid on
xlabel('År'); ylabel('Residual (M)')
legend('Linjär', 'Förfinad', 'Location','best')
title('Residualjämförelse')


% Delfråga c)

% Startgissning tas från del (b)
p = [d0; d1; d2; d3; k];              % p = [d0 d1 d2 d3 k]^T

maxit = 50;                           % max antal iterationer
tol   = 1e-10;                        % tolerans för parameterförändring

for it = 1:maxit
    d0 = p(1); d1 = p(2); d2 = p(3); d3 = p(4); k = p(5);

    % Modell och residual med nuvarande parametrar
    yhat3 = d0 + d1*t + d2*sin(k*t) + d3*cos(k*t);
    r3    = y - yhat3;

    % Jakobimatris J = ∂r/∂p för r = y - yhat
    J = [ -ones(n,1), ...
          -t, ...
          -sin(k*t), ...
          -cos(k*t), ...
          -t.*( d2.*cos(k*t) - d3.*sin(k*t) ) ];

    % Gauss–Newton-steget: lös J*delta ≈ -r (stabilt med '\')
    delta = J \ (-r3);

    % Enkel dämpning: backtracking om norm inte minskar
    s = 1;
    rnorm0 = norm(r3,2);
    while true
        p_try = p + s*delta;
        y_try = p_try(1) + p_try(2)*t + p_try(3)*sin(p_try(5)*t) + p_try(4)*cos(p_try(5)*t);
        if norm(y - y_try,2) <= rnorm0 || s < 1e-4
            break
        end
        s = 0.5*s;
    end

    p_new = p + s*delta;

    % Stoppvillkor
    if norm(p_new - p) <= tol*(norm(p)+eps)
        p = p_new; break
    end
    p = p_new;
end

% Slutparametrar och mått
d0g = p(1); d1g = p(2); d2g = p(3); d3g = p(4); kg = p(5);
Lg  = 2*pi/kg;                         % motsvarande period L
yhat3   = d0g + d1g*t + d2g*sin(kg*t) + d3g*cos(kg*t);
r3      = y - yhat3;
R2norm3 = norm(r3,2);

fprintf('\nGauss–Newton (k fri):\n');
fprintf('d0 = %.6f\nd1 = %.6f miljoner/år\n', d0g, d1g);
fprintf('d2 = %.6f\nd3 = %.6f\n', d2g, d3g);
fprintf('k  = %.6f 1/år  (L = %.2f år)\n', kg, Lg);
fprintf('||R_gn||_2 = %.6f\n', R2norm3);

% --- Sammanställ jämförelse ----------------------------------------------
figure
plot(1920+t, y, 'o'); hold on; grid on
plot(1920+t, yhat,  '-', 'LineWidth',1.5);        % linjär
plot(1920+t, yhat2, '-', 'LineWidth',1.5);        % förfinad med k från FFT
plot(1920+t, yhat3, '-', 'LineWidth',1.5);        % olinjär GN
xlabel('År'); ylabel('Miljoner invånare')
legend('Data','Linjär','Förfinad (k fix)','Olinjär (k fri)','Location','northwest')
title('Tre modeller: jämförelse')

figure
plot(1920+t, r,  '-o'); hold on; grid on
plot(1920+t, r2, '-');
plot(1920+t, r3, '-');
xlabel('År'); ylabel('Residual (M)')
legend('Linjär','Förfinad (k fix)','Olinjär (k fri)','Location','best')
title('Residualer för de tre modellerna')

% Utskrift: vilken norm är störst/minst och varför
labels = {'Linjär','Förfinad (k fix)','Olinjär (k fri)'};
vals   = [R2norm, R2norm2, R2norm3];
[~,imin] = min(vals); [~,imax] = max(vals);
fprintf('\nStörst ||R||_2: %s = %.6f\n', labels{imax}, vals(imax));
fprintf('Minst  ||R||_2: %s = %.6f\n', labels{imin}, vals(imin));
fprintf(['Orsak: fler frihetsgrader ger lägre minsta-kvadrat-residual; ', ...
         'modellen i (c) har flest (inkl. fritt k).\n']);

