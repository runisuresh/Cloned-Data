clc; clear all; close all;

%Generating the mechanical data from Ramachandra and Humphrey, J Biomech, 2019

%Define parameters
%Stress Free
Ro = 669D-6 / 2; %m
H = 66D-6; %m
Ri = Ro - H ; %m
L = 3.51D-3; %m
V = pi / 4 *((2*Ro)^2 - (2*Ri)^2) * L;
lz_iv = 1.51; %(-)
P_iv = 12.8*133.33; %in vivo pressure

%Material Parameters
%Eln, circ, axial, diag1, diag2
ce = 11.8953*1D3; %Pa
ct1 = 0.0751*1D3; %(-)
ct2 = 0.2805; %(-)
cz1 = 0.5828*1D3; %(-)
cz2 = 1.0321; %(-)
cd1 = 1.3762*1D3; %(-)
cd2 = 0.5358; %(-)
beta = [90 0 36.7467 360 - 36.7467] * pi/180;
mat.ce = ce;
mat.ck = [ct1 ct2 cz1 cz2 cd1 cd2 cd1 cd2];
mat.beta = beta;
mat.gk = [1 1 1 1];


P_range = (5:0.1:40) * 133.33;
nP = length(P_range);

lz_test = [0.95 * lz_iv  lz_iv 1.05 * lz_iv];
nZ = length(lz_test);

nVar = 7;
exp_data = zeros(nP,nVar,nZ);
for j = 1:nZ
    ro_range = zeros(nP,1);
    L_range = zeros(nP,1);
    lt_range = zeros(nP,1);
    sigmat_range = zeros(nP,1);
    sigmaz_range = zeros(nP,1);
    lz_range = zeros(nP,1);
    
    lz = lz_test(j);
    for i = 1:nP
        %Loads for the current step
        P = P_range(i);
        objective = @(ro)loaded_obj(ro, mat, P, Ri, H, V, L, lz);
        if i > 1
            ro_g = ro_range(i - 1);
        else
            ro_g = Ro;
        end 
        ro_range(i) = fzero(objective, ro_g);
        ri = sqrt(ro_range(i)^2 - V/(pi * lz * L));
        t = calc_stress(ro_range(i), mat, Ri, H, V, L, lz);
        lt_range(i) = (ri + ro_range(i)) / (2 *(Ri + H / 2));
        L_range(i) = calc_axial_force(ro_range(i), ri, t(3));
        sigmaz_range(i) = t(3);
        sigmat_range(i) = t(2);
        
        if (P == P_iv && j == 2)
            ri_iv = ri;
            h_iv = ro_range(i) - ri;
        end
        
    end
    
    exp_data(:,1,j) = P_range;
    exp_data(:,2,j) = ro_range;
    exp_data(:,3,j) = L_range;
    exp_data(:,4,j) = lt_range;
    exp_data(:,5,j) = lz;
    exp_data(:,6,j) = sigmat_range;
    exp_data(:,7,j) = sigmaz_range;
    
end

figure(1)
subplot(2,1,1)
hold on
plot(2 * exp_data(:,2,1)* 1D6, exp_data(:,1,1)/133, 'LineWidth', 1.5, 'Color',[0.75 0.75 0.75])
plot(2 * exp_data(:,2,2) * 1D6, exp_data(:,1,1)/133.33, 'LineWidth', 1.5, 'Color',[0.5 0.5 0.5])
plot(2 * exp_data(:,2,3) * 1D6, exp_data(:,1,1)/133.33, 'LineWidth', 1.5, 'Color',[0.25 0.25 0.25])
axis([600 1800 0 50])
xlabel('Outer Diameter \mum'); ylabel('Pressure (mmHg)')
set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold')
subplot(2,1,2)
hold on
plot(exp_data(:,4,1), exp_data(:,3,1)*1D3, 'LineWidth', 1.5, 'Color',[0.75 0.75 0.75])
plot(exp_data(:,4,2), exp_data(:,3,2)*1D3, 'LineWidth', 1.5, 'Color',[0.5 0.5 0.5])
plot(exp_data(:,4,3), exp_data(:,3,3)*1D3, 'LineWidth', 1.5, 'Color',[0.25 0.25 0.25])
% axis([1 1.8 0 30])
xlabel('Circ. Stretch'); ylabel('Force (mN)')
set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold')
    
figure(2)
subplot(2,1,1)
hold on
plot(exp_data(:,4,1), exp_data(:,6,1)/1D3, 'LineWidth', 1.5, 'Color',[0.75 0.75 0.75])
plot(exp_data(:,4,2), exp_data(:,6,2)/1D3, 'LineWidth', 1.5, 'Color',[0.5 0.5 0.5])
plot(exp_data(:,4,3), exp_data(:,6,3)/1D3, 'LineWidth', 1.5, 'Color',[0.25 0.25 0.25])
axis([0.5 2.5 0 160])
xlabel('Circ. Stretch'); ylabel('Circ. Stress (kPa)')
set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold')
subplot(2,1,2)
hold on
plot(exp_data(:,4,1), exp_data(:,7,1)/1D3, 'LineWidth', 1.5, 'Color',[0.75 0.75 0.75])
plot(exp_data(:,4,2), exp_data(:,7,2)/1D3, 'LineWidth', 1.5, 'Color',[0.5 0.5 0.5])
plot(exp_data(:,4,3), exp_data(:,7,3)/1D3, 'LineWidth', 1.5, 'Color',[0.25 0.25 0.25])
axis([0.5 2.5 0 140])
xlabel('Circ. Stretch'); ylabel('Axial Stress (kPa)')
set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'FontWeight', 'bold')

save('Ramachandra_LPA_exp_data.mat','exp_data', 'Ro', 'H', 'Ri', 'L', 'V', 'ri_iv', 'h_iv', 'lz_iv', 'P_iv');

%--------------------------------------------------------------------------
%Functions
function J = loaded_obj(ro, mat, P, Ri, H, V, L, lz)

    ri = sqrt(ro^2 - V/(pi * lz * L));
    h = ro - ri;
    lt = (ri + h / 2) / (Ri + H / 2);
    F = [1/(lt * lz); lt; lz];
    tex = constitutive(F, mat);
    
    lagrange = ones(3,1) * tex(1);
    t = tex - lagrange;

    t2_LP = P * ri / h;
    
    J = t(2) - t2_LP;

end

function tex = constitutive(F, mat)
    %Calculates the Caucy stress for a Fung-like material from the
    %Material properties, deformation, and strain
    
    %Elastin contributions
    S = ones(3,1) * mat.ce;
    
    %Fiber contributions
    for k = 1:4
        beta_k = mat.beta(k);
        lambda_k = sqrt((F(2)*sin(beta_k))^2 + (F(3)*cos(beta_k))^2);
        Q1 = (lambda_k^2 - 1);
        Q2 = Q1^2;
        c1 = mat.ck(2*(k - 1) + 1);
        c2 = mat.ck(2*(k - 1) + 2);
        g = mat.gk(k);
        
        for dir = 1:3
            
            switch dir
                case 1
                    G = 0;
                case 2
                    G = (g*sin(beta_k))^2;
                case 3
                    G = (g*cos(beta_k))^2;
            end
            
            S_k = c1 * Q1 * exp(c2 * Q2) * G;
            S(dir) = S(dir) + S_k;
        end
    end    
    
    tex = F.*S.*F;
end

function t = calc_stress(ro, mat, Ri, H, V, L, lz)

    ri = sqrt(ro^2 - V/(pi * lz * L));
    h = ro - ri;
    lt = (ri + h / 2) / (Ri + H / 2);
    F = [1/(lt * lz); lt; lz];
    tex = constitutive(F, mat);
    
    lagrange = ones(3,1) * tex(1);
    t = tex - lagrange;

end

function L = calc_axial_force(ro, ri, tz)

    L = pi*(ro^2 - ri^2) * tz;

end


