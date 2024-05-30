function tex = constitutive_mm4ff(F, mat)
    %microstructurally motivated 4 fiber family model
    %Calculates the Caucy stress for a Fung-like material from the
    %Material properties, deformation, and strain
    %weights by mass fraction for each constituent
    
    %Elastin contributions
    G = mat.ge.*mat.ge;
    S = ones(3,1) * mat.phie * mat.ce .* G;
    
    %Smooth muscle cell contribution
    phim = mat.phim;
    g = mat.gm;
    lambda_k = F(2) * g;
    if lambda_k < 1
        lambda_k = 1; 
    end
    Q1 = (lambda_k^2 - 1);
    Q2 = Q1^2;
    c1 = mat.cm(1);
    c2 = mat.cm(2);

    for dir = 1:3

        switch dir
            case 1
                G = 0;
            case 2
                G = g^2;
            case 3
                G = 0;
        end

        S_k = phim * c1 * Q1 * exp(c2 * Q2) * G;
        S(dir) = S(dir) + S_k;
    end
    
    %Fiber contributions
    phic = mat.phic;
    for k = 1:4
        beta_k = mat.beta(k);
        phi_k = mat.phik(k);
        g = mat.gc;
        lambda_k = g * sqrt((F(2) * sin(beta_k))^2 + (F(3) * cos(beta_k))^2);
        if lambda_k < 1
            lambda_k = 1; 
        end
        Q1 = (lambda_k^2 - 1);
        Q2 = Q1^2;
        c1 = mat.ck(2*(k - 1) + 1);
        c2 = mat.ck(2*(k - 1) + 2);
        
        for dir = 1:3
            
            switch dir
                case 1
                    G = 0;
                case 2
                    G = (g*sin(beta_k))^2;
                case 3
                    G = (g*cos(beta_k))^2;
            end
            
            S_k = phi_k * phic * c1 * Q1 * exp(c2 * Q2) * G ;
            S(dir) = S(dir) + S_k;
        end
    end    
    
    tex = F .* S .* F;
end

