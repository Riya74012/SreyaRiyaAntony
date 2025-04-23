clc;
clear;
syms x z w real

% Compute the db3 scaling function (phi) and wavelet function (psi)
[phi_vals, psi_vals, xval_grid] = wavefun('db3', 10);

% Define the scaling function phi(x) via interpolation
phi = @(xnum) interp1(xval_grid, phi_vals, xnum, 'linear', 0); 

hold on
fplot(phi,x)
hold off

% Settings
L = 5; 
r = 1;   
rho = L * r;  
A_mat = sym(zeros(rho, rho));  
X_values = [0.02447,0.2061,0.5,0.79389,0.97553]; 

% Construct the matrix A_mat
row_start = 1;  
for j = 1:L 
    block_rows = r;  
    block_cols = rho; 
    x_val = X_values(j);  
    A_block = sym(zeros(block_rows, block_cols));  
    
    for p = 0:block_rows-1
        for q = 0:block_cols-1
            term1 = phi(x_val - q);
            term2 = phi(x_val - q + rho);
            A_block(p+1, q+1) = term1 + term2 * z;
        end
    end
    
    A_mat(row_start:row_start + block_rows - 1, :) = A_block;
    row_start = row_start + r;  
end

disp('Symbolic Matrix A_mat:');
disp(A_mat);

% Calculate the determinant and inverse of A_mat
det_A = simplify(det(A_mat));
disp('Determinant of A_mat:');
disp(det_A);

if det_A ~= 0
    inv_A = simplify(inv(A_mat));
    disp('Inverse of A_mat:');
    disp(inv_A);

    % Compute symbolic Fourier transform of inverse matrix
    Fourier_inv_A = sym(zeros(size(inv_A)));
    [rA, cA] = size(inv_A);
    
    for i = 1:rA
        for j = 1:cA
            f_expr = subs(inv_A(i,j), z, exp(2*pi*1i*x));
            Fourier_inv_A(i,j) = int(f_expr * exp(-2*pi*1i * w * x), x, 0, 1);
        end
    end
    
    disp('Fourier transform of Inverse of A_mat:');
    disp(Fourier_inv_A);

    % Precompute Theta coefficients and shifts for numerical evaluation
    Theta_coeffs = cell(L, r);
    Theta_shifts = cell(L, r);

    for n = 0:L-1
        for p = 0:r-1
            coeffs = [];
            shifts = [];
            for q = 0:rho-1
                for nu = -20:20
                    coeff = limit(Fourier_inv_A(q+1, n*r + p + 1), w, nu);
                    if coeff ~= 0
                        coeffs(end+1) = coeff;
                        shifts(end+1) = rho*nu + q;
                    end
                end
            end
            Theta_coeffs{n+1, p+1} = coeffs;
            Theta_shifts{n+1, p+1} = shifts;
        end
    end

    % Plotting Theta functions numerically using phi()
    x_vals = linspace(-10, 10, 1000);
    colors = lines(L * r);
    figure;
    hold on;

    legend_entries = {};
    idx = 1;
    for n = 1:L
        for p = 1:r
            y_vals = zeros(size(x_vals));
            coeffs = Theta_coeffs{n, p};
            shifts = Theta_shifts{n, p};
            for k = 1:length(coeffs)
                y_vals = y_vals + double(coeffs(k)) * phi(x_vals - shifts(k));
            end
            plot(x_vals, y_vals, 'LineWidth', 2, 'Color', colors(idx, :));
            legend_entries{end+1} = sprintf('$\\Theta_{%d%d}$', n-1, p-1); %#ok<SAGROW>
            idx = idx + 1;
        end
    end

    xlabel('t');
    ylabel('$\Theta_{ni}(t)$', 'Interpreter', 'latex');
    %title('Plots of \Theta_{np}(x)', 'Interpreter', 'latex');
    legend(legend_entries, 'Interpreter', 'latex');
    grid on;
    hold off;

else
    disp('The matrix is singular and cannot be inverted.');
end
