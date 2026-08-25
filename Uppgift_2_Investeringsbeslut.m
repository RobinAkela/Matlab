% Parametrar
sigma  = 0.1;    
r      = 0.025; 
my     = 0.01;   
omega  = 1;     
kappa  = 4;     
lambda = 0.5;   
clear plot figure
set(0,'DefaultFigureVisible','on')


% Referensnivåer för σ -> 0 (används i verifikation av del (b))
wL = omega - r*lambda;      % w_L = ω - rℓ  (förväntad gräns för p_L)
wH = omega + r*kappa;       % w_H = ω + rκ   (förväntad gräns för p_H)

% Lös systemet för de givna parametrarna och skriv vias iterationsfel
% Kör Newton-solvern för givna parametrar. Funktionen returnerar även 'it' (antal steg),
% Men här bryr vi oss bara om pL och pH, därför tar vi bara de två första utgångarna.
[pL, pH] = solve_prices_sigma(sigma, my, r, omega, kappa, lambda, true);  % true => skriv ut ||x_{n+1}-x_n||
fprintf('Basfall σ=%.3f: pL=%.10f, pH=%.10f\n', sigma, pL, pH);
fprintf('Jämf. mot σ->0: pL-wL=%+ .3e, pH-wH=%+ .3e\n\n', pL - wL, pH - wH);

%  DEL (b): Verifiera pL->wL, pH->wH När σ Minskar
sigmas = [0.2 0.1 0.05 0.02 0.01 0.005];   % testvärden för σ (minskande)
pLvals = zeros(size(sigmas));              % vektor för p_L(σ)
pHvals = zeros(size(sigmas));              % vektor för p_H(σ)
fprintf('  sigma        pL            pH           pL-wL        pH-wH\n'); % Tabellhuvud

for i = 1:numel(sigmas) % Loop över σ
    s = sigmas(i);      % aktuellt σ
    [pLvals(i), pHvals(i)] = solve_prices_sigma(s, my, r, omega, kappa, lambda, false); % false => tyst
    fprintf('%7.4f   % .8f   % .8f   % .3e   % .3e\n', ...
        s, pLvals(i), pHvals(i), pLvals(i) - wL, pHvals(i) - wH);       % Tabellrad
end
fprintf('\n')

% (c) Svep i k i intervallet 0 till 5. Övriga parametrar fasta
K = linspace(0,5,201); % Intervallet mellan 0 till 5 med 201 punkter
pLc = zeros(size(K));                        % Resultatvektor
pHc = zeros(size(K));                        % Resultatvektor
for j = 1:numel(K)                           % Kör på intervallet K.
    kj = K(j);                               % Aktuellt K
    [pLc(j), pHc(j)] = solve_prices_sigma(sigma, my, r, omega, kj, lambda, false); % Löser för detta k (kj)
end

figure; plot(K, pLc, '-', K, pHc, '--'); grid on
xlabel('k'); ylabel('tröskelpris'); legend('p_L','p_H','Location','best')
title('Uppgift 2(c): p_L(k) och p_H(k)')

% (d) 1% osäkerhet i alla parametrar: experimentell störningsanalys
% Bas (k=4)
base = struct('sigma',sigma, 'my',my, 'r',r, 'omega',omega, 'k',kappa, 'l',lambda);
[pL0,pH0] = solve_prices_sigma(base.sigma, base.my, base.r, base.omega, base.k, base.l, false);
delta = 0.01;                          % Störningen på 1 %
names = {'sigma','my','r','omega','k','l'};
dL = zeros(6,1); 
dH = zeros(6,1);      % bidrag per parameter
for i = 1:6                              % högst 6 testkörningar + 1 bas = 7 körningar
    b = base;
    b.(names{i}) = b.(names{i})*(1+delta);
    [pL1,pH1] = solve_prices_sigma(b.sigma, b.my, b.r, b.omega, b.k, b.l, false);
    dL(i) = pL1 - pL0;                  % absolut förändring i pL
    dH(i) = pH1 - pH0;                  % absolut förändring i pH
end

% Sammanställ osäkerhet (Ger värsta fall och Root sum of squares (RSS) och rangordna bidrag
wcL = sum(abs(dL));   rssL = sqrt(sum(dL.^2));
wcH = sum(abs(dH));   rssH = sqrt(sum(dH.^2));

[~, idxH] = sort(abs(dH),'descend');    % mest till minst påverkan på pH
fprintf('\n(d) Bas: pL=%.10f, pH=%.10f\n', pL0, pH0)
fprintf('Osäkerhet pL: worst-case=%.3e, RSS=%.3e\n', wcL, rssL)
fprintf('Osäkerhet pH: worst-case=%.3e, RSS=%.3e\n\n', wcH, rssH)
fprintf('Bidrag per parameter (+1%%):\n')
fprintf('%8s   ΔpL         ΔpH        elasticitet_pL   elasticitet_pH\n', 'param')
for i=1:6
    theta0 = base.(names{i});
    eL = (dL(i)/pL0)/(delta);           % relativ känslighet (elasticitet)
    eH = (dH(i)/pH0)/(delta);
    fprintf('%8s  % .3e   % .3e    % .3e         % .3e\n', names{i}, dL(i), dH(i), eL, eH)
end
fprintf('\nRangordning störst→minst påverkan på pH: ')
fprintf('%s ', names{idxH}); 
fprintf('\n')

% Newton för givet σ 
function [pL, pH, it] = solve_prices_sigma(sigma, my, r, omega, kappa, lambda, printIter)
% [pL, pH, it] = ...   returnerar p_L och p_H samt antal iterationer 'it'.
% printIter=true  => skriv ut ||x_{n+1}-x_n|| varje iteration.

% Beräkna α<β som rötter till 0.5*σ^2 x^2 + (μ - 0.5*σ^2)x - r = 0
ab    = sort(roots([0.5*sigma^2, my - 0.5*sigma^2, -r]));  % sort => [α, β] i stigande ordning
alpha = ab(1);                                             % α
beta  = ab(2);                                             % β

% Bygg systemet F(x)=0 med x=[a; b; pL; pH]
% V0(p)=a p^β,  V1(p)=b p^α + p/(r-μ) - ω/r
F = @(x) [ ...
  x(1)*x(3)^beta - ( x(2)*x(3)^alpha + x(3)/(r-my) - omega/r ) - lambda;  % F1: V0(pL)-V1(pL)-ℓ
  x(1)*x(4)^beta - ( x(2)*x(4)^alpha + x(4)/(r-my) - omega/r ) + kappa;   % F2: V0(pH)-V1(pH)+κ
  x(1)*beta*x(3)^(beta-1) - ( x(2)*alpha*x(3)^(alpha-1) + 1/(r-my) );     % F3: V0'(pL)-V1'(pL)
  x(1)*beta*x(4)^(beta-1) - ( x(2)*alpha*x(4)^(alpha-1) + 1/(r-my) )      % F4: V0'(pH)-V1'(pH)
];

% Jacobimatris J = ∂F/∂(a,b,pL,pH) 
J = @(x) [ ...
  x(3)^beta, -x(3)^alpha, x(1)*beta*x(3)^(beta-1) - x(2)*alpha*x(3)^(alpha-1) - 1/(r-my), 0; % Rad 1
  x(4)^beta, -x(4)^alpha, 0, x(1)*beta*x(4)^(beta-1) - x(2)*alpha*x(4)^(alpha-1) - 1/(r-my); % Rad 2
  beta*x(3)^(beta-1), -alpha*x(3)^(alpha-1), x(1)*beta*(beta-1)*x(3)^(beta-2) - x(2)*alpha*(alpha-1)*x(3)^(alpha-2), 0; % Rad 3
  beta*x(4)^(beta-1), -alpha*x(4)^(alpha-1), 0, x(1)*beta*(beta-1)*x(4)^(beta-2) - x(2)*alpha*(alpha-1)*x(4)^(alpha-2)  % Rad 4
];

% Startgissning från σ=0-guiden (stabil och enkel)
wH = omega + r*kappa;                  % p_H vid σ=0
wL = omega - r*lambda;                 % p_L vid σ=0
x  = [1; 1; 0.9*wL; 1.1*wH];           % [a0; b0; pL0; pH0] nära referenserna
% Newtons metod: lös J(x)*dx = F(x), uppdatera x := x - dx
tau = 1e-10;                           % tolerans för ||x_{n+1}-x_n||_∞ (krav i (b))
Nmax = 100;                            % säkerhetsgräns
for it = 1:Nmax
    Fx = F(x);                         % residual F(x) i aktuell punkt
    dx = J(x) \ Fx;                    % lös linjärt system (dx är Newton-steget)
    x  = x - dx;                       % uppdatera till nästa punkt
    err = norm(dx, inf);               % ||x_{n+1}-x_n||_∞ = max(abs(dx))
    if printIter
        fprintf('%2d  ||x_{n+1}-x_n|| = %.3e\n', it, err); % uppgiftens utskrift
    end
    if err < tau, break; end           % stoppa när felet är < 1e-10 i alla komponenter
end

% Plocka ut lösningen (p_L, p_H) ur x
pL = x(3);                             % lägre tröskel
pH = x(4);                             % högre tröskel
end
