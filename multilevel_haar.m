function [V, Psi, P] = multilevel_haar(img,s);

    v = img(:);               % Vectorize the image
    N = length(img);          % Get the dimension of the image

    N2 = N^2;                 % Length of the image vector
    IN2 = speye(N2);          % A N^2 x N^2 identity matrix
    I2 = speye(2);            % A 2x2 identity matrix

    Psi = IN2;                % The complete wavelet transform as a matrix
    P = IN2;                  % The complete permutation matrix

    % Perform a multilevel wavelet decomposition by iterating over each level
    for n = 2.^(log2(N):-1:s+1)
        
        n2 = n^2;             % Length of the current sub-image vector
        I = speye(n/2);       % An n/2 x n/2 identity matrix
        
        % Generate a 1-D Haar wavelet transform matrix of size n x n
        A = 0.5*[ kron(I,[1,  1]); ...
                  kron(I,[1, -1])  ];

        % Convert the 1-D transform matrix to a 2-D transform of
        % size n^2 x n^2 which operates on the vectorized 2-D image
        Psi_n = kron(A,A);

        % Create a block matrix which allows us to apply the transform to only
        % the current sub-image, located at the top-left of the matrix
        Psi_n_N = IN2;
        Psi_n_N(1:n2,1:n2) = Psi_n;
        
        % Apply the transform to the image
        v = Psi_n_N*v;
        
        % Rearrange the transformed sub image to group the four coefficient
        % groups (ie Approximation, Horizontal, Vertical, and Diagonal) together
        P_n = kron(I2,kron([kron(I,[1 0]);kron(I,[0 1])],I));

        % Create a block matrix which allows us to apply the permutation to only
        % the current sub-image, located at the top-left of the matrix
        P_n_N = IN2;
        P_n_N(1:n2,1:n2) = P_n;
        
        % Apply the permutation to the image
        v = P_n_N*v;
        
        % Update Psi and P
        Psi = P_n_N * Psi_n_N * Psi;
        P = P_n_N * P;

    end
    
    % Permute the coefficient vector back to make it easier to visualize
    V = reshape(P'*v, N, []);
    
    % Remove the permutations from the Psi matrix
    Psi = P'*Psi;

end
