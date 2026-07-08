% Function to symmetrize the boundaries
function symm_array = symmetrize(array, pad_length)
    symm_array = [flipud(array(1:pad_length)); array; flipud(array(end-pad_length+1:end))];
end