function V = thresh(U,e)
    
    V = U;
    V(abs(U) < e) = 0;

end
