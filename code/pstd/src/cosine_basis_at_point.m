function phi = cosine_basis_at_point(x,L,N)

    phi = zeros(N,1);
    
    phi(1) = 1/sqrt(N);
    
    for m = 1:N-1
    
        phi(m+1) = ...
            sqrt(2/N) * cos(m*pi*x/L);
    
    end

end