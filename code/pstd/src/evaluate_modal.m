function u = evaluate_modal(Um,x,L)

    N = length(Um);
    
    u = Um(1)/sqrt(N);
    
    for m = 1:N-1
        u = u + ...
            sqrt(2/N) * Um(m+1) * cos(m*pi*x/L);
    end

end