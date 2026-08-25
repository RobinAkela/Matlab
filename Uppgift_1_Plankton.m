% Parametrar:
r = 0.3;
K = 108;
Rm = 0.7;
alpha = 5.7;
Z = 5;
tau = 1e-10;
Nmax = 1000;

% Funktioner och dess derivata
f = @(P) r.*P.*(1-(P./K)) - Rm.*Z.*((P.^2) ./ (alpha.^2+P.^2)); % f(p)
df = @(P) r - (2*r.*P)./K - Rm*Z .* (2*P.*alpha.^2) ./ (alpha.^2 + P.^2).^2; % f'(P)
g = @(P) ((P.^2)./K) + Rm*Z*(P.^2) ./ (r*(alpha^2+P.^2));  % g(p)
dg = @(P) 2*P./K + (Rm*Z./r) .* (2*P.*alpha.^2) ./ (alpha.^2 + P.^2).^2; % g'(p)


% a) Visar att (2) är fixpuntsinteration för f(P)=0
% Om härledningen är rätt gäller f(P) = r*(P - g(P)).
% Mäter max-avvikelse på ett tätt rutnär
Pgrid = linspace(1e-6,110,4000);
eqdiff = max(abs(f(Pgrid) - r.*(Pgrid - g(Pgrid)))); % Ska gå mot 0 
fprintf('största avvikelsen = %.3e\n', eqdiff)

% Hittar RÖTTER med fzero
% Tre bracket-intervalll där vi vet att det finns en rot i varje.
% fzero försöker hitta ett nollställe inom varje intervall.
I   = [3 5; 8 10; 90 100]; % bracket-intervall
rts = zeros(size(I,1),1);
for i = 1:size(I,1)
    rts(i) = fzero(f, I(i,:));
end
fprintf('Hittade rötter:%s\n', sprintf(' %.4f', rts)); % Utskrift av alla rötter

% Plotta f(P)
Pplot = linspace(0,110,2000);
figure; plot(Pplot, f(Pplot), 'LineWidth', 1);
grid on;
xlabel('P');
ylabel('f(P)');

% Gissningar efter undersökning av plottad graf:
x1 = 4;
x2 = 9.2;
x3 = 94;

% b) Teori för konvergens: Lokalt test
vals = abs(dg(rts)); % Beräknar absolutbeloppet av g' i varje rot
fprintf('|g''(r)| i rötter (i ordning) :%s\n', sprintf(' %.4f', vals));
for i = 1:numel(rts)
    status = 'instabil'; 
    if vals(i) < 1
        status = 'stabil';
    end
    fprintf('Rot %d: P = %.4f, |g''(r)| = %.4f -> %s\n', i, rts(i), vals(i), status);
end

% Välj en stabil rot att iterera mot (den första som uppfyller g' < 1)
idx = find(vals < 1, 1);
if isempty(idx)
    error('ingen stabil rot hittad') % Avbryter om alla rötter är instabila
end
stable_root = rts(idx); % De(n) stabila roten

% c) Själva fixpunktsiterationen
P = 0.9*rts(idx);                     % start nära men inte exakt i roten
err = nan(Nmax,1);
for n = 1:Nmax
    Pnext = g(P);                     % ett fixpunktssteg
    err(n) = abs(Pnext - P);          % felindikator
    disp([n, Pnext, err(n)])          % tabell: [steg, P_n, |P_{n+1}-P_n|]
    if err(n) < tau, break; end       % stoppa på tolerans
    P = Pnext;
end
P_star = Pnext;
steps = n;
fprintf('Fixpunkt ~ %.12g, iterationer = %d\n', P_star, steps);

% d) Kvottest
err = err(1:steps);
if numel(err) >= 2
    rq = err(2:end)./err(1:end-1);    % err_(k+1})/err_k
    fprintf('sista kvot = %.4f, |g''(P*)| = %.4f\n', rq(end), ...
        abs(dg(P_star)));
end

% NEWTONS METOD (e–g)

% Startgissningar
starts = [x1, x2, x3];

% Behållare för Newton-resultat
newton_roots = zeros(numel(starts),1); % rot per start
newton_steps = zeros(numel(starts),1); % steg per start

% e) Kör Newtons metod för alla tre rötter
for j = 1:numel(starts)
    Pn = starts(j); % Startvärde för denna rot
    fprintf('\nNewton från start %.4f (rot %d):\n', Pn, j);
    fprintf('n      P_n             |P_{n+1}-P_n|\n');
    errN = nan(Nmax,1); % Lagra stegfel för delfråga f)

    for k = 1:Nmax
        % Newton-steg: P_(n+1) = f(P_n) / f'(P_n)
        step = f(Pn) ./ df(Pn);      % Steglängd f/df
        Pn1 = Pn - step;             % Nästa punkt
        errN(k) = abs(Pn1 - Pn);     % Felindikator
        fprintf('%-5d  %.12g    %.3e\n', k, Pn1, errN(k));

        if errN(k) < tau % Stoppa på samma tolerans
            break
        end
        Pn = Pn1; % Fortsätt iterering
    end
newton_roots(j) = Pn1;  % Hittad rot
newton_steps(j) = k;    % Antal steg

% f) Tumregelkontroll (kvadratisk)
% Om Newton är kvadratisk: err_(k+1) ungefär lika med C*err_k^2
% Vi visar detta via sista kvoten q = err_(k+1)/err_k^2
errN = errN(1:k);
if numel(errN) >= 2
    q_last = errN(end) / (errN(end-1)^2); % Bör vara ungefär konstant
    digits_prev = -log10(errN(end-1));    % Ungefär antal korrekta siffor före
    digits_last = -log10(errN(end));      % Och efter
    fprintf('kvot err_{k+1}/err_k^2 ≈ %.3e, siffror: %.2f -> %.2f\n', ...
        q_last, digits_prev, digits_last);
end
end

% Sammanfattning för alla tre Newton-körningar
fprintf('\nNewtons rötter:%s\n', sprintf(' %.4f', newton_roots));
fprintf('Newtons steg för ovan i samma ordning  :%s\n', sprintf(' %d', ...
    newton_steps));

% g) Jämför fixpunkt vs Newton kring stabila roten
% Samma start i båda: 0.9 * stabil rot från fixpunktdelen
start_compare = 0.9*stable_root;

% Fixpunkt: Kör snabb miniloop enbart för stegräkning
P = start_compare; steps_fix = 0;
for k = 1:Nmax
    Pn1= g(P);
    steps_fix = k;
    if abs(Pn1 - P) < tau
        break
    end
    P = Pn1;
end

% Newton: kör från samma start
P = start_compare; steps_newt = 0;
for k = 1:Nmax
    Pn1 = P - f(P)/df(P);
    steps_newt = k;
    if abs(Pn1 - P) < tau
        break
    end
    P = Pn1;
end

% Utskrift av jämförelsen
fprintf('\nJämförelse från start %.4f kring stabil rot:\n', start_compare);
fprintf('Fixpunkt steg: %d (linjär, faktor ≈ %.3f)\n', steps_fix, ...
    abs(dg(stable_root)));
fprintf('Newton   steg: %d (kvadratisk)\n', steps_newt);

% Svara på delfråga g)
if steps_newt < steps_fix
    fprintf(['Slutsats (g): Newton är snabbast eftersom konvergensen är kvadratisk, ' ...
        'fixpunkt är linjär (faktor ≈ %.3f).\n'], ...
        abs(dg(stable_root))); 
    % abs(dg(stable_root) visar hur mycket felet minskar per steg nära
    % fixpunkten och förklarar varför fixpunkt är långsammare än Newton.
else
    fprintf(['Slutsats (g): Fixpunkt tog färre steg här, ' ...
        '        men metoden är linjär (faktor ≈ %.3f) medan Newton är kvadratisk.\n'], ...
        abs(dg(stable_root))); 
end






