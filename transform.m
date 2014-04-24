function V = transform(H,U)

    V = reshape(H*U(:),length(U),[]);

end
