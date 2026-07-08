function data = init_FDTD(len_x, c, dt, dh, space_order, alpha1, alpha2)

    laplacian = get_laplacian_kernel(space_order);

    data.laplacian = laplacian;
    data.c = c;
    data.dt = dt;
    data.dh = dh;
    data.alpha1 = alpha1;
    data.alpha2 = alpha2;

end