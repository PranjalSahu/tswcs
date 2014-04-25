function Phi = sampling_matrix(M,N,s)

    Phi = zeros(M,N);
    Phi(1:s,1:s) = eye(s);
    
    for n = s+1:N
        Phi(s+1:M,n) = round(1.2*rand(M-s,1)+0.4)-1;
        
    end
    % Phi = zeros(M,N);
    % Phi(1:s,1:s) = eye(s);
    
    % for n = s+1:N
    %     Phi(s+1:M,n) = 2*floor(2*rand(M-s,1))-1;
    % end

end
