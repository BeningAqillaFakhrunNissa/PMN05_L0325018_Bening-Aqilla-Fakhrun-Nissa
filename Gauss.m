%Bening Aqilla Fakhrun Nissa - L0325018
%Eliminasi Gauss
%Soal: 2x1 + x2 - x3 = 3 ; 4x1 + 3x2 + x3 = 9 ; -2x1 + x2 + 2x3 = 4
%   A = [2 1 -1; 4 3 1; -2 1 2];
%   b = [3; 9; 4];
%   x = Gauss(A, b)

function x = Gauss(A, b)
tic
[n, l] = size(A);
fprintf('Matriks augmented awal:\n');
disp([A b]);
for i = 1 : n-1
    if A(i, i) == 0
        error('Pivot nol pada baris %d', i);
    end
    for h = i+1 : n
        m = A(h, i) / A(i, i);
        fprintf('m%d%d = %g\n', h, i, m);
        A(h, :) = A(h, :) - m * A(i, :);
        b(h, :) = b(h, :) - m * b(i, :);
    end
end
fprintf('Hasil forward elimination:\n');
disp([A b]);
x(n, :) = b(n, :) / A(n, n);
for i = n-1 : -1 : 1
    x(i, :) = (b(i, :) - A(i, i+1:n) * x(i+1:n, :)) / A(i, i);
end
fprintf('Hasil substitusi mundur:\n');
for i = 1 : n
    fprintf('x%d = %.6f\n', i, x(i));
end
toc
end
