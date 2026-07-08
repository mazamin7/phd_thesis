function simulation_parameters = get_simulation_parameters()
        msg_method_left = "Choose the update method for left";
        opts_method = ["FDTD" "Fourier"];
        choice_method_left = menu(msg_method_left, opts_method);
    
        msg_method_right = "Choose the update method for right";
        choice_method_right = menu(msg_method_right, opts_method);

        msg_space_order = "Choose the spatial FDTD/residual order";
        opts_space = ["2" "4" "6" "8"];
        choice_space = menu(msg_space_order, opts_space);

    space_order = 2 * choice_space;

    simulation_parameters.method_left = choice_method_left;
    simulation_parameters.method_right = choice_method_right;
    simulation_parameters.space_order = space_order;
end
