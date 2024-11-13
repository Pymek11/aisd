%% 1
A=[3,1;
   2,1;
   1,1];
b=[3;
    4;
    5];
x= inv(A' * A)*(A' * b),

%%  2

close all;
clear all;

A=[3,1;
   2,1;
   1,1];
b=[3;
    4;
    5];
[Q,R]=qr(A),
Qt_b =Q' *b,

n = size(A,2),
r2 = Qt_b(1:n),

R1=R(1:n,:);
x = R1 \ r2,


%% 3
% approx_line.m
% Regresja liniowa: y = a*x + b
clear all; close all;
% W wyniku pomiaru otrzymano nastepujace liczby ( x = numer pomiaru, y = wartosc )
x = [ 1 2 3 4 5 6 10 ];
y = [ 0.912 1.013 1.035 1.057 1.062 1.082 1.097 ];
figure; plot( x, y, 'b*' ); title('y=f(x)'); grid; pause
% Aproksymacja linia prosta: y = a * x + b
if(0) % OGOLNIE - rozwiazanie rownania macierzowego
xt = x'; yt = y'; N = length(xt); % X * ab
X = [ xt, ones(N,1) ]; % y(n) = a*x(n) + b = | x(1) 1 | * | a |
ab = X \ yt; % y(1) = a*x(1) + b = | x(2) 1 | | b |
a = ab(1), b = ab(2), % y(2) = a*x(2) + b = | x(3) 1 |
else % W TYM PRZYPADKU - na podstawie wyprowadzonych wzorow
xm = mean( x ); % srednia wartosc wektora x
ym = mean( y ); % srednia wartosc wektora y
xr = x - xm; % wektor x - srednia x (od kazdego elementu)
yr = y - ym; % wektor y - srednia y (od kazdego elementu)
a = (xr * yr') / (xr * xr') % obliczenie wsp a prostej, to samo
% inaczej: a = sum( xr .* yr ) / sum( xr .* xr )

b = ym - a * xm % obliczenie wsp b prostej
end
%figure; plot( x, y, 'b*', x, a*x+b, 'k-' ); title('y=f(x)'); grid; pause
% Takze wielomiany wyzszych rzedow
p = polyfit( x, y, 4 ),  a=p(1), b=p(2), c= p(3), 
figure; plot( x, y, 'b*', x, polyval(p,x), 'r-' ); title('y=f(x)'); grid; pause

%% 4

% approx_channel.m
% Estymacja odpowiedzi impulsowej kanalu transmisji, u nas karty dzwiekowej
clear all; close all;
h = [3; -2; 1 ]; % symulowana odpowiedz impulsowa kanalu
%load h.dat % rzeczywista odpowiedz impul
% sowa kanalu
L = length(h); % liczba wag
SNR_values = [0,5,10,15];
N_values = [30, 40, 50];
for snr = SNR_values
    for N = N_values
        prbs = 2*round( rand(N,1) )-1; % pobudzenie kanalu
        %prbs=[-1,-1,-1,-1,1,-1,1,-1,1,1,1,-1,1,1,-1,1,-1,1,-1,-1,1]’; % krotkie
        %prbs=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]’; % test dla X
        x = [ prbs ];

        r = x(L:-1:1); % pierwszy wiersz
        c = x(L:N); % pierwsza kolumna
        X = toeplitz(c,r), pause % macierz pobudzenia o strukturze Toeplitza
        y = X*h; % wynik przejscia pobudzenia przez uklad h(n)
        y_ = awgn(y,snr); % dodanie szumu, tak aby otrzymac zadany SNR
        %he = pinv(X)*y; % czyli inv(X’*X) * X’*y = Rxx * rxy
        he = X \ y_; % to samo zoptymalizowane obliczeniowo
        plot(1:L,h,'ro',1:L,he,'bx'); title('Zadane i obliczone h(n)'); 
        title(['SNR = ', num2str(snr),'dB, Sygnal = ',num2str(N)]);
        legend('Rzeczywista wartosc', 'Estymowana wartosc')
        grid;
    end
end

%% 5
% approx_krata.m
close all; clear all;
% Generacja/wczyranie obrazka
N = 512; Nstep = 32;
[img, cmap] = imread('Lena512.bmp'); img = double(img); % Lena
%img = zeros(N,N); % czarny kwadrat
if(0) % opcjonalna biala siatka
for i=Nstep:Nstep:N-Nstep, img(i-1:i+1,1:N) = 255*ones(3,N); end
for j=Nstep:Nstep:N-Nstep, img(1:N,j-1:j+1) = 255*ones(N,3); end
end
imshow( img, cmap ); pause
% Dodawanie znieksztalcen beczkowych
a = [ 1.06, -0.0002, 0.000005 ]; % wspolczynniki wielomianu znieksztalcen
x=1:N; y=1:N; 
cx = N/2+0.5; % SRODEK DEFORAMACJI
cy = N/2+0.5;
[X,Y] = meshgrid( x, y ); % wszystkie x,y
r = sqrt( (X-cx).^2 + (Y-cy).^2 ); % wszystkie odleglosci od srodka
R = a(1)*r.^1 + a(2)*r.^2 + a(3)*r.^3; % zmiana odleglosci od srodka
Rn = R ./ r; % normowanie
imgR = interp2( img, (X-cx).*Rn+cx, (Y-cy).*Rn+cy ); % interploacja
figure;
subplot(1,2,1),imshow(img,cmap); title('Oryginal');
subplot(1,2,2),imshow(imgR, cmap); title('Rybie oko'); pause
% Estymacja znieksztalcen beczkowych
i = Nstep : Nstep : N-Nstep; j=i; % polozenie linii w pionie i poziomie
[I,J] = meshgrid( i, j ); % wszystkie (x,y) punktow przeciec
r = sqrt( (I-cx).^2 + (J-cy).^2 ); % wszystkie promienie od srodka
R = a(1)*r + a(2)*r.^2 + a(3)*r.^3; % odpowiadajace punkty obrazu znieksztalconego
r = sort( r(:) ); % sortowanie
R = sort( R(:) ); % sortowanie
aest1 = pinv([ r.^1, r.^2, r.^3 ])*R; aest1 = [ aest1(end:-1:1); 0]; % rozw.1
aest2 = polyfit( r, R, 3)'; % rozw.2
[ aest1, aest2 ], pause % porownanie
aest = aest1; % wybor rozwiazania
% Wielomian R=f(r) i odwrotny r=g(R)
r = 0:N/2; % wybrane promienie
R = polyval( aest, r); % R=f(r) wielomianu znieksztalcen
figure; subplot(121); plot(r,R), title('R=f(r)');
ainv = polyfit( R, r, 3), % wspolczynniki wielomiany odwrotnego
subplot(122); plot(R,r), title('r=g(R)'); pause
% Korekta znieksztalcen beczkowych
[X,Y] = meshgrid( x, y ); % wszystkie punkty (x,y) znieksztalconego
R = sqrt( (X-cx).^2 + (Y-cy).^2 ); % wszystkie zle promienie
Rr = polyval( ainv, R ); % wszystkie dobre promienie
Rn = Rr./R; % normowanie
imgRR=interp2( imgR, (X-cx).*Rn+cx, (Y-cy).*Rn+cy ); % interpolacja
figure;
subplot(1,2,1),imshow(imgR,cmap); title('Wejscie - efekt rybie oko');
subplot(1,2,2),imshow(imgRR,cmap); title('Wyjscie - po korekcie');
colormap gray
pause
