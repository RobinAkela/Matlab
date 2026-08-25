function u4_a_build_matrix
a = 1; b = 2;
Tv = 400; Tl = 22;
n = 4;

[A,f,r,h] = build_A_f(n,a,b,Tv,Tl);

fprintf('h = %.5f\n',h);
disp('r-inre:'); disp(r.');
disp('A ='); disp(A);
disp('f ='); disp(f);
end

%% hjälpfunktion
function [A,f,rInterior,h] = build_A_f(n,a,b,Tvatska,Tluft)
h = (b-a)/(n+1);
rInterior = a + (1:n)'*h;

A = zeros(n); f = zeros(n,1);

for k=1:n
    rk = rInterior(k);
    alpha =  rk/h^2 - 1/(2*h);
    beta  = -2*rk/h^2;
    gamma =  rk/h^2 + 1/(2*h);

    if k>1
        A(k,k-1) = alpha;
    else
        f(k) = f(k) - alpha*Tvatska;
    end

    A(k,k) = beta;

    if k<n
        A(k,k+1) = gamma;
    else
        f(k) = f(k) - gamma*Tluft;
    end
end
end
