% Laddar eiffel1
% (Fx,Fy) = (+1,-1) alltså nedåt höger med beloppet sqrt(2)

load eiffel1.mat
n = numel(xnod); % Antalet noder n, behöver för att dimensionera b.

% Välj nod "j" och bygg b så att kraften är:
j = 1; 
b = zeros(2*n,1);
b(2*j-1) = 1; % Fx_j
b(2*j) = -1;  % Fy_j

x = A\b;                  % Lös systemet: Ax = b med backslash. Ger förskjutningar delta-x och y
xbel = xnod + x(1:2:end); % Nya x-koordinater x^bel = x+delta(x)
ybel = ynod + x(2:2:end); % Nya y-koordinater x^bel = y+delta(y)

% Figur för plottning av torn
figure;
hold on;
axis equal;
grid on;

trussplot(xnod,ynod,bars);                  % original-truss
trussplot(xbel,ybel,bars);                  % belastad-truss
plot(xnod(j),ynod(j),'k*','MarkerSize',9);  % markera vald nod med asterix

% Mät tid för Ax=b för var och en av de sex modellerna
% Använd godtyckligt b(randn). Använder även tic/toc. 
% Plotta tid mot N=2n i loglog
% Jämför därefter visuellt med teori ~ O(N^3) för tät Gausseliminering.

files = {'eiffel1.mat';'eiffel2.mat';'eiffel3.mat';'eiffel4.mat';'eiffel5.mat';'eiffel6.mat'};

% Snabbtest av kod
% Om alla filer ska köras, sätt:
% files_b = files; 
% files_c = files;
files_b = {'eiffel1.mat';'eiffel2.mat';'eiffel3.mat';'eiffel4.mat';'eiffel5.mat';'eiffel6.mat'};
files_c = {'eiffel1.mat';'eiffel2.mat';'eiffel3.mat'};


Nvals = zeros(numel(files_b),1);    % Lagra N = 2*n per modell för x-axeln
tSolve = zeros(numel(files_b),1);   % Lagra uppmätt tid för Ax=b per modell för y-axeln
rng(0)                              % slump för b, ska vara godtycklig.

for k = 1:numel(files_b)                          % gå igenom alla modeller
    S = load(files_b{k},'A','xnod');              % ladda A och xnod
    Ak = S.A;  nk = numel(S.xnod);  Nk = 2*nk;    % antal obekanta
    Nvals(k) = Nk;                                % spara N
    rhs = randn(Nk,1);                            % godtyckligt b

    tic                                           % start tid
    xk = Ak\rhs;                                  % backslash (Gauss) – resultatet ej vidare använt
    tSolve(k) = toc;                              % stanna tid
end


figure;                            % Ny figur för tidsdiagrammet
loglog(Nvals,tSolve,'o-');         % loglog-plot
grid on;
hold on;
xlabel('Antal obekanta  N = 2n')   % Kopplar x-axeln till definitionen i texten.
ylabel('Tid för Ax=b  (s)')        % Visar vad som mätts.
title('Uppgift 3(b): backslash-tid per modell')


N0 = Nvals(1);                                                    % Väljer första punkten som ref. för skalning
t0 = tSolve(1) + eps;                                             % Ref. tid, +eps undviker 0 vid små tider
loglog(Nvals, t0*(Nvals/N0).^3, '--')                             % Rita ~N^3. Visar den teoretiska lutningen för tät eliminering.
legend('Mätt tid','Ref ~ N^3','Location','northwest')             % Gör jämförelsen läsbar i figuren

% Sammanställningstabell med kolumner för modell, N ,j*, T* och tider för
% fyra metoderna.
Res = table(strings(numel(files_c),1), zeros(numel(files_c),1), zeros(numel(files_c),1), zeros(numel(files_c),1), ...
            zeros(numel(files_c),1), zeros(numel(files_c),1), zeros(numel(files_c),1), zeros(numel(files_c),1), ...
            'VariableNames',{'Modell','N','j_star','T_star','Naiv_s','LU_s','Gles_s','GlesLU_s'});

for k = 1:numel(files_c)                                               % Gå igenom alla eiffel modeller.
    S = load(files_c{k},'A','xnod','ynod','bars');                     % Ladda A och noddata
    A = S.A; xnod = S.xnod; ynod = S.ynod; 
    n = numel(xnod); N = 2*n;                                        % antal obekanta
    Res.Modell(k) = string(files_c{k}); 
    Res.N(k) = N;
    
    if k==1                                  % Visa gles/bandad struktur (ej i tid)
        figure; spy(A); title('spy(A) – kontroll av gles/bandad');
    end

% (i) Naiv backslash i loopen
    T = zeros(n,1);                           % T_j = ||x_j|| ska beräknas för alla j.
    t0 = tic;                                 % Tidmätning för beräkningarna.
    for j = 1:n                               % Loopa över alla noder j
        b = zeros(N,1); b(2*j-1)=1; b(2*j)=-1;% Bygg b_j med Fx=1, Fy=-1 i nod j.
        x = A\b;                              % lös Ax=b_j med backslash.
        T(j) = norm(x);                       % Definiera T_j som normen av x_j
    end
    t_naiv = toc(t0);                         % Total tid för naiv metod över alla j
    [Tstar,jstar] = max(T);                   % Käsnligaste noden j* är den som maximerar T_j

% Rimlighetsplot för eiffel1 (inte i tidsmätning)
    if k==1
    figure; 
    axis equal;
    grid on;
    hold on
    trussplot(xnod,ynod,S.bars);
    plot(xnod(jstar),ynod(jstar),'r*','MarkerSize',8);
    title(sprintf('Känsligaste nod (naiv), %s: j^*=%d', files_c{k}, jstar));
    legend('Truss','Känsligaste nod');
    end

% (ii) Full LU: faktor en gång + lös alla b_j
t0 = tic; [L,U,P] = lu(A); t_fact_full = toc(t0); % Gör LU-faktorisering en gång; tiden ska ingå.
Tlu = zeros(n,1);
t0 = tic;                                         % Mät lösningstiden separat.
for j = 1:n
    b = zeros(N,1); b(2*j-1)=1; b(2*j)=-1;        % Samma b_j som ovan
    x = U\(L\(P*b));                              % Full LU
    Tlu(j) = norm(x);
end
t_LU = t_fact_full + toc(t0);                     % Total tid för full LU
[~, j_lu] = max(Tlu);
    
% (iii) Gles backslash
As = sparse(A);                                     % Konvertera till gles, ska ej tidsmätas
Tsp = zeros(n,1);                                   
t0 = tic;                                           % Mät endast lösningarna
for j = 1:n
    b = zeros(N,1); b(2*j-1)=1; b(2*j)=-1;          % Samma lastdefinition
    x = As\b;                                       % Backslash på gles matris
    Tsp(j) = norm(x);                               % T_j = ||x_j||, Används för att hitta känsligaste noden.

end
t_gles = toc(t0);                                   % Total tid för gles slashback
[~, j_sp] = max(Tsp);

% (iv) Gles + LU
t0 = tic; [Ls,Us,Ps,Qs] = lu(As); t_fact_sp = toc(t0); % Gles LU, P*A*Q = L*U faktoriseingstid ska ingå.
Tglu = zeros(n,1);
t0 = tic;                                              % mät lösningstiden separat.

for j = 1:n
    b = zeros(N,1); b(2*j-1)=1; b(2*j)=-1;  % Samma b_j
    y = Ls\(Ps*b); y2 = Us\y; x = Qs*y2;    % Fram- och baksubstitution samt permuteringsåterställning.
    Tglu(j) = norm(x); 
end 
t_glesLU = t_fact_sp + toc(t0);             % Total tid för gles LU = faktor + lösningar.
[~, j_glu] = max(Tglu);

% Verifiering av:
% Samma känsligaste nod och storlek på maximal förskjutning. 
assert(jstar==j_lu && jstar==j_sp && jstar==j_glu, 'Olika känsligaste nod.');
assert(max(abs([Tstar-max(Tlu), Tstar-max(Tsp), Tstar-max(Tglu)])) < 1e-9, 'Olika T*.');

% Spara rad
Res.Modell(k) = string(files_c{k}); % fyll tabellens rad för aktuell modell
Res.N(k)      = N;                  % antal obekanta
Res.j_star(k) = jstar;              % känsligaste nodens index
Res.T_star(k) = Tstar;              % motsvarande euklidisk norm T* 
Res.Naiv_s(k) = t_naiv;             % tid för naiv metod
Res.LU_s(k)   = t_LU;               % tid för full LU
Res.Gles_s(k) = t_gles;             % tid för gles backslasj
Res.GlesLU_s(k)= t_glesLU;          % tid för gles + LU
end

disp(Res)

% Varför är LU snabbare? 
% Problemet: lös Ax = b_j för många b_j

% Naiv kör A\b för varje b_j
% Medan LU faktoreiserar A = LU en gång och återanvänder
% Sedan dessa för två snabbare steg: L*y = P*b_j och U*x_j = y
% Tät matris har en tidskostnad på ungefär:
% Naiv metod: ≈ (2/3)*N^3
% LU ≈ engångskostnad på ((2/3)*N^3) sedan 2*N^2 per b
% Gles per b ≈ nnz(L) + nnz(U) där nnz = number o nonzeros i matrisen
% b ≈ nnz(L) + nnz(U) är ofta mindre än N^2.
% Slutsats: Att göra den dyra faktoriseringen en gång och återanvända den
% gör LU betydligt snabbare när du har många b_j.

% Vilken metod löser snabbast?

% För enstaka högerled väljer A\b en bra metod utoamtiskt; för glesa A gör
% Matlab gles LU i bakgrunden.
% För många högerled är Gles + LU snabbast, Full LU näst snabbast
% Därefter gles backslash i loop och sist naiv backslash i loop.

% För vilken modell blir tidsvinsten störst? Varför?

% De största modellerna (med störst N) och de mest glesa/bandade
% Det vill säga Naiv metod som växer med ungefär N^3 för varje lösning.
% Lu med gleshet gör en tung faktorisering en gång och får därefter snabba
% lösningar, totaltiden växer långsammare.

