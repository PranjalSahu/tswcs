function V = my_inverse_haar(img,s);

V = img;
N = length(V);

for n = 2.^(s+1:log2(N))
    A = [ kron(eye(n/2),[1,  1]); ...
          kron(eye(n/2),[1, -1])  ];
    V(1:n,1:n) = A'*V(1:n,1:n)*A;
end

end
