function t = calc_stress(ro, mat, Ri, H, V, L, lz, F_iv)

    ri = sqrt(ro^2 - V/(pi * lz * L));
    h = ro - ri;
    lt = (ri + h / 2) / (Ri + H / 2);
    F = [1/(lt * lz); lt; lz] ./ F_iv;
    tex = constitutive_mm4ff(F, mat);
    
    lact = 1; %lt / F_iv(2);
    C = mat.CB - mat.CS * 0;
    parab_act = 1 - ((mat.lm - lact) / (mat.lm - mat.l0))^2;
    act_stress = mat.phim * mat.T_act * (1 - exp(-C^2)) * lact * parab_act;
    %display( act_stress )
    
    lagrange = ones(3,1) * tex(1);
    t = tex - lagrange;
    t(2) = t(2) + act_stress;

end

