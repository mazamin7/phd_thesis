function laplacian = get_laplacian_kernel(space_order)

    assert(space_order == 2 || space_order == 4 || space_order == 6 || space_order == 8);
    n = space_order/2;

    if n == 1
        laplacian = [1, -2, 1];
    elseif n == 2
        laplacian = [-1/12, 4/3, -5/2, 4/3, -1/12];
    elseif n == 3
        laplacian = [1/90, -3/20, 3/2, -49/18, 3/2, -3/20, 1/90];
    else
        laplacian = [-1/560, 8/315, -1/5, 8/5, -205/72, 8/5, -1/5, 8/315, -1/560];
    end

end
