%% 1.1
clear all; close all;
m=5; v0=50; alpha=30; h=50; g=9.81;
alpha = alpha/180*pi;
x = 0 : 1 : 350;
y = h + tan(alpha)*x - g / (2*v0^2*cos(alpha)) * x.^2;
figure; plot(x,y); xlabel('x'); ylabel('y'); title('y(x)'); grid;
hold on;
zasieg = fzero(@(x) h + tan(alpha) * x - g / (2 * v0^2 * cos(alpha)) * x.^2, [0, 350]),
plot(zasieg, 0, 'rp', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % funkcja fzero bierze y=0 wstawia i oblicza x
                                                                 % plottowanie punktu y=0         

                                                                 

%% 1.2 z oporem powietrza dla kuli armatniej 
clear all; close all;
m=5; v0=50; alpha=30; h=50; g=9.81; 
C=0.47; % wspolczynnik powietrza\
p=1.225;% gestosc popwietrza
r= 0.1;% zalozone 10cm promien kuli ( wzor wzuety z internetu)(wartosci powietrza tez)
A= pi * r ^2;
b = @(v) 1/2 * C * p * A* v;% b zmienic na wspolczynnik oporu powietrza
alpha = alpha/180*pi;
x = 0 : 1 : 350;
y = (tan(alpha) + m * g ./ (b(v0) * v0 * cos(alpha))) .* x + g * (m^2) ./ (b(v0)^2) .* log(1 - x .* b(v0) ./ (m * v0 * cos(alpha)));
figure; plot(x,y); xlabel('x'); ylabel('y'); title('y(x)'); grid;
hold on;
zasieg = fzero(@(x) (tan(alpha) + m * g ./ (b(v0) * v0 * cos(alpha))) .* x +g * (m^2) ./ (b(v0)^2) .* log(1 - x .* b(v0) ./ (m * v0 * cos(alpha))), [10, 350]),
plot(zasieg, 0, 'rp', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % funkcja fzero bierze y=0 wstawia i oblicza x
                                                                 % plottowanie punktu y=0         
%% 2

clear all; close all;
it = 15;
a = pi - pi/5; b = pi + pi/5; 
f = @(x) sin(x); 
fp = @(x) cos(x); 


cb = nonlinsolvers(f, fp, a, b, 'bisection', it);
tolerance = pi * 0.00001; % Dokładność 0.001%
diff = abs(cb - pi);
iter_needed = find(diff < tolerance, 1);

% Wyświetlenie wyników
fprintf('Potrzebne iteracje: %d\n', iter_needed);
disp('Końcowe oszacowania:');
disp(cb(1:iter_needed));

%% 2.12

clear all; close all;
it = 10;
a_par = -1; b_par = 1; c_par = 0; % dla  45 a = -1 b=3 c=-2 (x0 = 1)
                                  % dla  5  a = -1 b =2.087 c =-1.087
                                  % dla  80 a = -1 b =7.67 c=-6.67
f = @(x) a_par*x.^2 + b_par*x + c_par;%% funckja
fp = @(x) 2*a_par*x + b_par;% pochodna funckji

cb = nonlinsolvers(f, fp, -2, 2, 'bisection', it);

disp('Końcowe oszacowania dla paraboli:');
disp(cb);

%% 2.21
it = 12;
a = pi-pi/5; b=pi+pi/5; 
f = @(x) sin(x); % funkcja
fp = @(x) cos(x); % pochodna funkcji
x = 0 : 0.01 : 2*pi;
plot( x, f(x), 'b-', x, fp(x),'r-'); grid; xlabel('x'); title('f(x), fp(x)');
legend('Funkcja','Jej pochodna'); pause
cb = nonlinsolvers( f, fp, a, b, 'bisection', it );
cr = nonlinsolvers( f, fp, a, b, 'regula-falsi', it);
cn = nonlinsolvers( f, fp, a, b, 'newton-raphson', it);
plot( 1:it,cb,'o-', 1:it,cr,'*', 1:it,cn,'^-'); xlabel('iter'); title('c(iter)')
grid on, legend('Bisection','Regula-Falsi','Newton-Raphson');


%% 3.1
% Test dla metody regula-falsi i funkcji sinus
clear all; close all;
it = 10; 
a = pi - pi/5; b = pi + pi/5; % Zakres początkowy
f = @(x) sin(x); % Funkcja
fp = @(x) cos(x); % Pochodna funkcji


cr = nonlinsolvers(f, fp, a, b, 'regula-falsi', it);

tolerance = pi * 0.00001; % Dokładność 0.001%(tolerancja bledu)
diff = abs(cr - pi);
iter_needed = find(diff < tolerance, 1);
fprintf('Regula-falsi: Potrzebne iteracje: %d\n', iter_needed);
disp('Końcowe oszacowania:');
disp(cr(1:iter_needed));
%% Test dla paraboli i metody regula-falsi
clear all; close all;
it = 50;
a = -2; b = 2;
a_par = -1; b_par = 3; c_par = -2; % dla  45 a = -1 b=3 c=-2 (x0 = 1)
                                   % dla  5  a = -1 b =2.087 c =-1.087
                                   % dla  80 a = -1 b =7.67 c=-6.67
f = @(x) a_par*x.^2 + b_par*x + c_par;
fp = @(x) 2*a_par*x + b_par;

% Test metody regula-falsi
cr = nonlinsolvers(f, fp, a, b, 'regula-falsi', it);

% Wyświetlenie wyników
disp('Końcowe oszacowania dla paraboli:');
disp(cr);

%% 9.5

% Test dla metody newton-raphson i funkcji sinus
clear all; close all;
it = 10; 
a = pi - pi/5; b = pi + pi/5; % Zakres początkowy
f = @(x) sin(x); % Funkcja
fp = @(x) cos(x); % Pochodna funkcji


cr = nonlinsolvers(f, fp, a, b, 'newton-raphson', it);

tolerance = pi * 0.00001; % Dokładność 0.001%(tolerancja bledu)
diff = abs(cr - pi);
iter_needed = find(diff < tolerance, 1);
fprintf('metoda newton-raphson: Potrzebne iteracje: %d\n', iter_needed);
disp('Końcowe oszacowania:');
disp(cr(1:iter_needed));

%% 9.5.2
% Test dla paraboli i metody newton-raphson
clear all; close all;
it = 50;
a = -2; b = 2;
a_par = -1; b_par = 3; c_par = -2; % dla  45 a = -1 b=3 c=-2 (x0 = 1)
                                   % dla  5  a = -1 b =2.087 c =-1.087
                                   % dla  80 a = -1 b =7.67 c=-6.67
f = @(x) a_par*x.^2 + b_par*x + c_par;
fp = @(x) 2*a_par*x + b_par;

% Test metody newtona - raphsona
cr = nonlinsolvers(f, fp, a, b, 'newton-raphson', it);

% Wyświetlenie wyników
disp('Końcowe oszacowania dla paraboli:');
disp(cr);


%% 9.6
% equnonlin_newtonraphson.m
clear all; close all; format long;
A = 2;
f = @(x) (x - A)^(1+x); % wzor funkcji
fp = @(x) 2*x; % wzor jej pochodnej
x(1) = 10 ;
x(2) = 9;
K=10;
% punkty startowe (estymata poczatkowa), liczba iteracji
for k = 3 : K
    x(k)=x(k-1) - f(x(k-1)) / fp(x(k-1)); % z wzorem pochodnej
    % x(k)=x(k-1) - f(x(k-1))/((f(x(k-1))-f(x(k-2)))/(x(k-1)-x(k-2))); % z estym. pochodnej
end
blad = sqrt(A) - x(K),
plot(1:K, x,'b.-', 1:K, sqrt(A)*ones(1,K),'r-'); grid; title('x(k)'); pause

