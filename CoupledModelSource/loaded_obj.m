function J = loaded_obj(ro, mat, P, Ri, H, V, L, lz, F_iv)

    t = calc_stress(ro, mat, Ri, H, V, L, lz, F_iv);
    ri = sqrt(ro^2 - V/(pi * lz * L));
    h = ro - ri;
    t2_LP = P * ri / h;
    
    J = t(2) - t2_LP;

end


