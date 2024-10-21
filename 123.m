%1
% Definicja wektorów
u = [1; 2; 3];
v = [4; 5; 6];

% Definicja macierzy
A = [1, 2, 3; 4, 5, 6; 7, 8, 9];
B = eye(3);

% 1. Dodawanie wektorów
u_plus_v = u + v;

% 2. Mnożenie wektora przez liczbę skalarową
scalar_mult_u = 2 * u;

% 3. Iloczyn skalarny wektorów
dot_product_uv = dot(u, v);

% 4. Mnożenie macierzy przez wektor
A_times_u = A * u;

% 5. Transpozycja wektora
u_transposed = u';

% 6. Mnożenie macierzy A i B
A_times_B = A * B;

% 7. Macierz odwrotna do B (która jest jednostkowa, więc odwrotna to sama B)
inv_B = inv(B);

% 8. Dodanie rzędu do macierzy A
new_row = [10, 11, 12];
A_extended = [A; new_row];

% 9. Wyznacznik macierzy A
det_A = det(A);

% 10. Macierz odwrotna do A (jeśli istnieje)
inv_A = NaN;  % Sprawdzanie, czy macierz jest odwracalna
if det_A ~= 0
	inv_A = inv(A);
else
	disp('Macierz A nie jest odwracalna');
end

% 11. Macierz przeciwna do A
neg_A = -A;

% Wyświetlanie wyników
disp('u + v = '), disp(u_plus_v);
disp('2 * u = '), disp(scalar_mult_u);
disp('u . v = '), disp(dot_product_uv);
disp('A * u = '), disp(A_times_u);
disp('u transposed = '), disp(u_transposed);
disp('A * B = '), disp(A_times_B);
disp('inv(B) = '), disp(inv_B);
disp('Macierz A z dodanym rzędem: '), disp(A_extended);
disp('Wyznacznik A = '), disp(det_A);
if ~isnan(inv_A)
	disp('Macierz odwrotna do A = '), disp(inv_A);
end
disp('Macierz przeciwna do A = '), disp(neg_A);

%2 
% matrix_obwod_dc.m
clear all; close all;
R1 = 10; R2 = 20; R3 = 30; R0 = 40;
E1 = 1; E2 = 2; E3 = 3;
A = [ R1+R2, -R2, 0; ...
-R2, R2+R3, -R3; ...
0, -R3, R3+R0 ],
b = [ E1-E2; ...
E2-E3; ...
E3 ],
% x = ?
x1 = inv(A)*b; % inv(A) = A^(-1)
x2 = pinv(A)*b; % pinv(A) = (A^T * A)^(-1) * A^T
x3 = A \ b; % minimaliacja bledu sredniokwadratowego
% Metoda Cramera
for k=1:length(b)
Ak = A; Ak(:,k) = b; % (w,k) = (:,k)
x4(k) = det( Ak ) / det(A);
end
x4 = x4.’;
[ x1, x2, x3, x4 ], pause
%b - natężenie 

%3 
% Definicja częstotliwości i elementów obwodu
f = 50; % Hz
Ze = 10 + j*5; % Impedancja Ze
Zs1 = 20; % Rezystor Zs1
Zs2 = 20; % Rezystor Zs2
Cd = 1e-6; % Pojemność kondensatora w faradach
Zd = 1 / (j*2*pi*f*Cd); % Impedancja kondensatora Zd
Lo = 0.1; % Indukcyjność w henrach
Zo = j*2*pi*f*Lo; % Impedancja induktora Zo
E = 230; % Napięcie źródła E

% Tworzenie macierzy równań
A = [ (1/Ze + 1/Zs1), -1/Zs1, 0; ...
  	-1/Zs1, (1/Zs1 + 1/Zs2 + 1/Zd), -1/Zs2; ...
  	0, -1/Zs2, (1/Zs2 + 1/Zo) ];

b = [ E/Ze; 0; 0 ]; % tylko jedno źródło napięcie więc jedna wartość E/Z

% Rozwiązywanie równań potencjałów węzłowych
U = A \ b;

% Wyświetlanie wyników
U1 = U(1);
U2 = U(2);
U3 = U(3);

disp('Potencjały węzłów:');
disp(['U1 = ', num2str(U1), ' V']);
disp(['U2 = ', num2str(U2), ' V']);
disp(['U3 = ', num2str(U3), ' V']);

