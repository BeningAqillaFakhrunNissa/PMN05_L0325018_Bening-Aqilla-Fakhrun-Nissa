%Bening Aqilla Fakhrun Nissa - L0325018
%Dekomposisi LU
%Soal: 2x1 + x2 - x3 = 3 ; 4x1 + 3x2 + x3 = 9 ; -2x1 + x2 + 2x3 = 4
%   A = [2 1 -1; 4 3 1; -2 1 2];
%   b = [3; 9; 4];
%   x = dekomposisiLU(A, b)
function x = dekomposisiLU(A, b)
tic;
fprintf('Hasil dekomposisi LU\n')
[L, U, b] = LU_fwdelim(A, b);
fprintf('Matriks L:\n'); disp(L);
fprintf('Matriks U:\n'); disp(U);
fprintf('Pembuktian L*U (harus sama dengan A):\n'); disp(L * U);
fprintf('Selisih A - L*U:\n'); disp(A - L * U);
[n, m] = size(L);
z = zeros(n, 1);
x = zeros(n, 1);
z(1) = b(1) / L(1, 1);
for i = 2 : n
    z(i) = (b(i) - L(i, 1:i-1) * z(1:i-1)) / L(i, i);
end
fprintf('Hasil Ly = b (forward substitution), y =\n'); disp(z);
x(n) = z(n) / U(n, n);
for i = n-1 : -1 : 1
    x(i) = (z(i) - U(i, i+1:n) * x(i+1:n)) / U(i, i);
end
fprintf('Hasil Ux = y (backward substitution):\n');
for i = 1 : n
    fprintf('x%d = %.6f\n', i, x(i));
end
toc;
end

%Forward elimination untuk Dekomposisi LU (menyimpan pengali m ke L)
function [L, U, b] = LU_fwdelim(A, b)
[n, l] = size(A);
L = zeros(n);
for i = 1 : n-1
    if A(i, i) == 0
        error('Pivot nol pada baris %d', i);
    end
    for h = i+1 : n
        m = A(h, i) / A(i, i);
        L(h, i) = m;
        A(h, :) = A(h, :) - m * A(i, :);
    end
end
L = L + eye(n);
U = A;
end
