clear all; close all;

load test_images
U = double(test_image{4})/256-0.5;
u = U(:);
N = length(u);

[V, Psi, P] = multilevel_haar(U,0);

display('Generating Sampling Matrix...');
Phi = sampling_matrix(2000,N);

H = Phi*Psi';

display('Random Sampling...');
y = H*u;

display('Basis Pursuit...')
cvx_begin
    variable t(N)
    minimize( norm(t,1) )
    subject to
        y == Phi*t
cvx_end

display('Complete!')

T = reshape(t,size(U));
showme(T)
