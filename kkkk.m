%% zadanie 1.1

close all;clear;

% Dane dla wielomianu f(x) = ax + b
a = 2; % współczynnik przy x
b = 3; % wyraz wolny
f = @(x) a*x + b; % wielomian

% Punkty i krok
x = [0, 3]; % punkty
h = x(2) - x(1); % krok

% Dokładna całka i pochodna
calka_matlab = integral(f, x(1), x(2)); % dokładna całka w przedziale [x(1), x(end)]
pochodna = a; % dokładna pochodna f'(x)

% 1. Całka metodą trapezów
calka_wzor = (x(2) - x(1)) / 2 * (f(x(1)) + f(x(2)));

% 2. Pochodna metodą różnic wprzód
pochodna_wzor = (f(x(2)) - f(x(1))) / h;

% Wyświetlenie wyników
disp(['Dokładna całka: ', num2str(calka_matlab)]);
disp(['Całka metodą trapezów: ', num2str(calka_wzor)]);
disp(['Dokładna pochodna: ', num2str(pochodna)]);
disp(['Pochodna metodą różnic wprzód: ', num2str(pochodna_wzor)]);


%% zadanie 1.2

close all; clear;

% Dane
x = [0, pi/8, pi/4]; % punkty x
y = sin(x); % wartości funkcji sin(x)

% 1. Całka metodą trapezów na przedziale [0, pi/4]
h_c = x(end) - x(1); % szerokość przedziału
trapez = (h_c / 2) * (y(1) + y(end));
disp(['Całka metodą trapezów: ', num2str(trapez)]);

% 2. Iloraz różnicowy wprzód dla x=0, h=pi/8
h_r= x(2) - x(1); % krok
iloraz = (y(2) - y(1)) / h_r;
disp(['Pochodna metodą ilorazu różnicowego wprzód: ', num2str(iloraz)]);

% sin(x) można przybliżyć w okolicy zera za pomocą szeregu taylora, również
% można użyć interpolacji lagrangea a puzniej aproksymacji

%% zadanie 2.1

close all; clear;

% f sinus
xz = [0,pi/4,pi/2];

syms x;
y=sin(x);

yz = [sin(xz(1)),sin(xz(2)),sin(xz(3))];

pochodna=cos(x);

%dokladne wartosci pochodnych z sinx
p = [cos(xz(1));cos(xz(2));cos(xz(3))];
h=pi/4;

px1=(1/(2*h))*(-3*yz(1)+4*yz(2)-yz(3));
px2=(1/(2*h))*(yz(3)-yz(1));
px3=(1/(2*h))*(yz(1)-4*yz(2)+3*yz(3));

err0=abs(px1-p(1));
err1=abs(px2-p(2));
err2=abs(px3-p(3)); % otrzymane blezy sa duże

%% zadanie 2.2

close all; clear;

%f wielomianowa
xz= [1,2,3];

syms x;

y=0.5+x+2*(x^2);

yz = [0.5+xz(1)+2*(xz(1)^2),0.5+xz(2)+2*(xz(2)^2),0.5+xz(3)+2*(xz(3)^2)];

pochodnay_y=4*x+1; % pochodna z funkcji wielomianowej

pz= [4*xz(1)+1,4*xz(2)+1,4*xz(3)+1];

% obliczenie pochodnej na podstawie 3 wzorow
h1=1;
ppx1=(1/(2*h1))*(-3*yz(1)+4*yz(2)-yz(3));
ppx2=(1/(2*h1))*(yz(3)-yz(1));
ppx3=(1/(2*h1))*(yz(1)-4*yz(2)+3*yz(3));

%obliczenie wartosci błędów
erx0=abs(ppx1-pz(1));
erx1=abs(ppx2-pz(2));
erx2=abs(ppx3-pz(3)); % błędy wynoszą 0

% wielomian  y = 0.5+x +2x^2+3x^3

wielomian=0.5+x+2*(x^2)+3*(x^3); %skorzystamy z wczesniejszego xz

ww = [0.5+xz(1)+2*(xz(1)^2)+3*(xz(1)^3),0.5+xz(2)+2*(xz(2)^2)+3*(xz(2)^3),0.5+xz(3)+2*(xz(3)^2)+3*(xz(3)^3)];

pochodna=9*(x^2)+4*x+1;

wp = [9*(xz(1)^2)+4*xz(1)+1,9*(xz(2)^2)+4*xz(2)+1,9*(xz(3)^2)+4*xz(3)+1];

%obliczenie pochodnej na podstawie 3 wzorow
pw1=(1/(2*h1))*(-3*ww(1)+4*ww(2)-ww(3));
pw2=(1/(2*h1))*(ww(3)-ww(1));
pw3=(1/(2*h1))*(ww(1)-4*ww(2)+3*ww(3));

ex1=abs(pw1-wp(1));
ex2=abs(pw2-wp(2));
ex3=abs(pw3-wp(3)); %bledy mniejsze niz 10, duze

%% Zadanie 3

clear; close all;

N = 1000;
x = sin(2*pi/N*(0:N-1));

% Szum o sile 
strength = [0.01, 0.1, 0.5];
err_d1 = zeros(length(strength), N);
err_d2 = zeros(length(strength), N);

for s = 1:length(strength)
    x_n = x + strength(s) * randn(1, N);
    
    % Pochodne numeryczne
    % Pochodne numeryczne
    for k = 1:N
        if k == 1  % Różnica przednia dla pierwszego punktu
            d1_numerical = (x_n(k+1) - x_n(k));                   % Różnica przednia
            d2_numerical = (x_n(k+1) - 2*x_n(k) + x_n(N)) / 1^2;  % Przednia różnica dla 2. pochodnej
        elseif k == N  % Różnica tylna dla ostatniego punktu
            d1_numerical = (x_n(k) - x_n(k-1));                   % Różnica tylna
            d2_numerical = (x_n(1) - 2*x_n(k) + x_n(k-1)) / 1^2;  % Tylna różnica dla 2. pochodnej
        else  % Różnice centralne dla pozostałych punktów
            d1_numerical = (x_n(k+1) - x_n(k-1)) / 2;               % Różnica centralna dla 1. pochodnej
            d2_numerical = (x_n(k+1) - 2*x_n(k) + x_n(k-1)) / 1^2;  % Centralna różnica dla 2. pochodnej
        end
        
        % Pochodne rzeczywiste
        d1_actual = cos(2*pi/N*((k+1)));
        d2_actual = -sin(2*pi/N*((k+1)));
        
        % Średnie błędy przybliżenia
        err_d1(s, k) = mean(abs(d1_numerical - d1_actual));
        err_d2(s, k) = mean(abs(d2_numerical - d2_actual));
    end
end

% najmnuejszy blad
[min_error_d1, idx_d1] = min(err_d1, [], 2);
[min_error_d2, idx_d2] = min(err_d2, [], 2);

disp('Wyniki dla pierwszej pochodnej:');
disp(['Najmniejszy błąd dla różnych wartości szumu: ', num2str(min_error_d1')]);
disp(['Odpowiadające wartości k: ', num2str(idx_d1')]);

disp('Wyniki dla drugiej pochodnej:');
disp(['Najmniejszy błąd dla różnych wartości szumu: ', num2str(min_error_d2')]);
disp(['Odpowiadające wartości k: ', num2str(idx_d2')]);

%% zadania 4

clear; close all;

% Parametry
N = 1000; % Liczba próbek
T = 1; % Okres sinusa
fs = N/T; % Częstotliwość próbkowania
t = (0:N-1)/fs; % Czas
x = sin(2*pi*t); % Oryginalny sygnał (sinus)
strengt = 0.1; % Siła szumu

% Zaszumiony sygnał
x_n = x + strengt * randn(1, N);

% Obliczanie pierwszej pochodnej numerycznie (oryginalny sygnał)
d1_x = zeros(1, N);
for k = 2:N-1
    d1_x(k) = (x(k+1) - x(k-1))/2; % Różnice centralne
end
d1_x(1) = (x(2) - x(1)); % Różnica przednia dla pierwszego punktu
d1_x(N) = (x(N) - x(N-1)); % Różnica tylna dla ostatniego punktu

% Obliczanie pierwszej pochodnej numerycznie (zaszumiony sygnał)
d1_x_n = zeros(1, N);
for k = 2:N-1
    d1_x_n(k) = (x_n(k+1) - x_n(k-1))/2; % Różnice centralne
end
d1_x_n(1) = (x_n(2) - x_n(1)); % Różnica przednia dla pierwszego punktu
d1_x_n(N) = (x_n(N) - x_n(N-1)); % Różnica tylna dla ostatniego punktu

% Filtracja FIR (odszumianie sygnału)
% Filtr FIR, liczba wag 7, pasmo przepustowe 0.2 (5/25)
h = fir1(7, 5/20); % Filtr FIR
x_n_filtered = filter(h, 1, x_n); % Odszumiony sygnał

% Obliczanie pierwszej pochodnej dla odszumionego sygnału
d1_x_n_filtered = zeros(1, N);
for k = 2:N-1
    d1_x_n_filtered(k) = (x_n_filtered(k+1) - x_n_filtered(k-1))/2; % Różnice centralne
end
d1_x_n_filtered(1) = (x_n_filtered(2) - x_n_filtered(1)); % Różnica przednia
d1_x_n_filtered(N) = (x_n_filtered(N) - x_n_filtered(N-1)); % Różnica tylna

% Wykres porównania pochodnych
figure;
plot(t, d1_x, 'b', 'LineWidth', 1.5); hold on;
plot(t, d1_x_n, 'r', 'LineWidth', 1.5);
plot(t, d1_x_n_filtered, 'g', 'LineWidth', 1.5);
legend('Pochodna oryginalnego sygnału', 'Pochodna zaszumionego sygnału', 'Pochodna odszumionego sygnału');
xlabel('Czas (s)');
ylabel('Pochodna');
title('Porównanie numerycznych pierwszych pochodnych');
grid on;

%% zadanie 5

clear; close all;

K = 28; % Nowa liczba
w = 0 : pi/100 : pi;

d1 = 1/12 * [-1 8 0 -8 1];
d2 = firls(K-1, [0 0.5 0.7 1], [0 0.5*pi 0 0], 'differentiator');
d3 = firpm(K-1, [0 0.5 0.7 1], [0 0.5*pi 0 0], 'differentiator');

figure;
plot(w/pi, abs(freqz(d1, 1, w))/pi, 'b-', ...
     w/pi, abs(freqz(d2, 1, w))/pi, 'r--', ...
     w/pi, abs(freqz(d3, 1, w))/pi, 'm-.');
xlabel('f/fnorm'); title('|D(fnorm)|'); grid;
legend('DIFF', 'LS', 'MIN-MAX');

% Zastosowanie wag filtrów do sygnału sinusoidalnego
f = 0.1; 
t = 0:0.01:1;
x = sin(2*pi*f*t);

y1 = filter(d1, 1, x); % filtracja różniczkująca skończona
y2 = filter(d2, 1, x); % filtracja za pomocą firls
y3 = filter(d3, 1, x); % filtracja za pomocą firpm

figure;
%subplot(3,1,1); plot(t, x, 'b-', t, y1, 'r--'); title('Filtracja DIFF');
%subplot(3,1,2); plot(t, x, 'b-', t, y2, 'r--'); title('Filtracja LS');
%subplot(3,1,3); plot(t, x, 'b-', t, y3, 'r--'); title('Filtracja MIN-MAX');
