function Phi = sampling_matrix(M,N)

    Phi = zeros(M,N);

    for n = 1:N
        Phi(:,n) = 2*floor(2*rand(M,1))-1;
    end

end
