function u3_a_euler_springs
% Uppgift 3(a): Systemform + testkörning med Framåt Euler.

eta  = 0.05;
beta = 0.10;

u0 = [-1; 0; -3; 0];   % [x; x'; y; y']
t0 = 0;
T  = 20;
h  = 0.01;             % valfri fin steglängd

[t,U] = euler_system(@f_springs, t0, T, u0, h, eta, beta);

x = U(1,:); y = U(3,:);

figure;
plot(t,x,'b',t,y,'r'); grid on;
xlabel('t'); ylabel('position');
legend('x(t)','y(t)');
title('Fjädersystem – Framåt Euler (uppgift 3a)');

fprintf('Euler: x(20)=%.4f, y(20)=%.4f (h=%.4f)\n', x(end), y(end), h);
end

%%  HJÄLPFUNKTIONER 
function dudt = f_springs(t,u,eta,beta)
% u = [x; x'; y; y'] = [x; v; y; w]
x = u(1); v = u(2); y = u(3); w = u(4);
dxdt = v;
dvdt = -eta*v - 2*x - beta*x^3 + y - 1;
dydt = w;
dwdt = -eta*w + 0.5*(x - y) - 1;
dudt = [dxdt; dvdt; dydt; dwdt];
end

function [t,U] = euler_system(f,t0,T,u0,h,eta,beta)
N = round((T - t0)/h);
t = linspace(t0,T,N+1);
U = zeros(length(u0),N+1);
U(:,1) = u0;

for k = 1:N
    U(:,k+1) = U(:,k) + h * f(t(k), U(:,k), eta, beta);
end
end
