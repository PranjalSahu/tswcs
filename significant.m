function V = significant(U,e)

    V = U;
    V(abs(U) < e) = 0;
    V(find(V)) = 1;

end
