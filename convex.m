clear all;

load test_images

SQRTN = length(test_image{1});
N = SQRTN*SQRTN;

display('Generating Transform Matrix...');
[~, Psi, ~] = multilevel_haar(speye(SQRTN),0);
    
for M = [ 2000 4000 6000 ]

    display('Generating Sampling Matrix...');
    Phi = sampling_matrix(2000,N,0);
    H = Phi*Psi';

    for i = 1:length(test_image)
        U = double(test_image{i})/256-0.5;
        u = U(:);
        N = length(u);
        V = transform(Psi.',U);

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
        UU = reshape(transform(Psi,T),size(U));

        fig1 = figure(1); clf;
        title1 = sprintf('Image Estimate, %i Samples, Basis Pursuit',M);
        subplot(2,2,2); showme(UU); title(title1); colorbar off;
        title2 = sprintf('True Image: "%s"',image_title{i});
        subplot(2,2,1); showme(U); title(title2); colorbar off;
        title3 = sprintf('Coefficient Estimate');
        subplot(2,2,4); showme(significant(T,0.5)); title(title3); colorbar off;
        title4 = sprintf('True Wavelet Coefficients');
        subplot(2,2,3); showme(significant(V,0.5)); title(title4); colorbar off;
        filename = sprintf('out/%s_m%i_cvx.eps',image_title{i},M);
        print(fig1, '-depsc2', filename);
    end
end
