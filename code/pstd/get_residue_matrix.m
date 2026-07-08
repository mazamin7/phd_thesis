function C = get_residue_matrix(N_x, space_order)

    assert(N_x == 2*floor(N_x/2), 'N_x is odd');

    assert(space_order == 2 || space_order == 4 || space_order == 6 || space_order == 8);
    n = space_order/2;

    if n == 1
        coefs = [1];
    elseif n == 2
        coefs = [-1/12, 4/3];
    elseif n == 3
        coefs = [1/90, -3/20, 3/2];
    else
        coefs = [-1/560, 8/315, -1/5, 8/5];
    end

    C = sparse(N_x,N_x);

    for i = 1:n
        row = N_x/2-n+i;
        
        col1 = N_x/2-i+1;
        col2 = N_x/2+i;

        vals = [-coefs(1:i) coefs(i:-1:1)];
        C(row,col1:col2) = vals;

        row = N_x/2+n-i+1;

        C(row,N_x/2-i+1:N_x/2+i) = -vals;
    end

end
