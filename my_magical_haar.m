function V = my_magical_haar(img);

V = img;
N = length(V);

for n = 2.^(log2(N):-1:1)
    A = [ kron(eye(n/2),[1,  1]); ...
          kron(eye(n/2),[1, -1])  ];
    V(1:n,1:n) = 0.25*A*V(1:n,1:n)*A';
end

end
