function HN = myhaar(N);

H2 = [1 1; 1 -1]/sqrt(2);

HN = H2;
for i = 2:log2(N)
    HN = [kron(HN,[1,1]); kron(eye(2^(i-1)),[1,-1])]/sqrt(2);
end

end
