% viscosity modified law taken from Pries et al Resistance to blood flow in
% microvessel in vivo.
% P is terminal branch pressure
function hemo = calculate_hemo_morphometric_tree_termBC(Q, P)

plasma_viscosity=0.0124; % Poise

global thickness

load morphometric_tree.mat

 n_orders = max(seg_size(:,3));
 number_seg=size(seg_connectivity,1);    

 WSS_tree=zeros(number_seg,1);
 Sigma_tree=zeros(number_seg,1);
 Q_tree=zeros(number_seg,1);
 r_tree = zeros(number_seg,1);
 P_tree = zeros(number_seg,2); %upstream and downstream pressure of each segment
 dP_tree = zeros(number_seg,1);
 n_seg_ord = zeros(n_orders,1);
 P_seg_ord = zeros(n_orders,2);
 Q_seg_ord = zeros(n_orders,2);
 Sigma_seg_ord = zeros(n_orders,2);
 WSS_seg_ord = zeros(n_orders,2);
 
 Q_tree(1)=Q;

 d=seg_size(1,1)*2*10000; %convert to micron
 viscosity=rel_viscosity(d)*plasma_viscosity;
 WSS_tree(1)=4*viscosity*Q/(pi*seg_size(1,1)^3);
 r_tree(1) = 8*viscosity*seg_size(1,2)/(pi*(seg_size(1,1))^4);
 dP_tree(1) = Q_tree(1) * r_tree(1);
 
 for j=1:number_seg

     if seg_connectivity(j,2)>0 && seg_connectivity(j,3)>0
                
         %indices for left and right child branches
         ind_l = seg_connectivity(j,2);
         ind_r = seg_connectivity(j,3);
         
         %add to current count of segments
         seg_order_l = seg_size(ind_l,3);
         seg_order_r = seg_size(ind_r,3);
         n_seg_ord(seg_order_l) = n_seg_ord(seg_order_l) + 1;
         n_seg_ord(seg_order_r) = n_seg_ord(seg_order_r) + 1;

         %Store resistances for calculating flow split
         rl=seg_input_resistance(ind_l);
         rr=seg_input_resistance(ind_r);

         %Store geometry for local resistance calc
         radius_l=seg_size(ind_l,1);
         radius_r=seg_size(ind_r,1);

         length_l=seg_size(ind_l,2);
         length_r=seg_size(ind_r,2);

         %Update flows
         %Note: Uses full downstream resistance from other child
         Q_tree(ind_l)=Q_tree(j)*rr/(rl+rr);
         Q_tree(ind_r)=Q_tree(j)*rl/(rl+rr);

         %Find WSS and downstream pressure for left child
         d=radius_l*2*10000;%convert to micron
         viscosity=rel_viscosity(d)*plasma_viscosity;
         WSS_tree(ind_l)=4*viscosity*Q_tree(ind_l)/(pi*(radius_l)^3);
         r_tree(ind_l) = 8*viscosity*length_l/(pi*(radius_l)^4);
         dP_tree(ind_l) = Q_tree(ind_l) * r_tree(ind_l);

         %Find WSS and downstream pressure for right child
         d=radius_r*2*10000;%convert to micron
         viscosity=rel_viscosity(d)*plasma_viscosity;
         WSS_tree(ind_r)=4*viscosity*Q_tree(ind_r)/(pi*(radius_r)^3);
         r_tree(ind_r) = 8*viscosity*length_r/(pi*(radius_r)^4);
         dP_tree(ind_r) = Q_tree(ind_r) * r_tree(ind_r);

     end

 end
 
 P_term = P * 133.33 * 10; %dynes / cm^2
 P_tree(number_seg, 2) = P_term;

 %Go from terminal to top to find root pressure
 seg = number_seg;
 while (seg > 1)
    par_seg = seg_connectivity(seg,1);

    %Assign upstream pressure from daughter to downstream pressure
    %of parent
    P_tree(seg, 1) = P_tree(seg, 2) + dP_tree(seg);
    P_tree(par_seg, 2) = P_tree(seg, 1);
    
    seg = par_seg;
     
 end
 P_tree(seg, 1) = P_tree(seg, 2) + dP_tree(seg);
 
 %Go down the tree again and find pressures using root pressure
 for j=1:number_seg

     if seg_connectivity(j,2)>0 && seg_connectivity(j,3)>0
         
         %indices for left and right child branches
         ind_l = seg_connectivity(j,2);
         ind_r = seg_connectivity(j,3);
         
         %Assign downstream pressure from parent to upstream pressure
         %of child
         P_tree(ind_l, 1) = P_tree(j, 2);
         P_tree(ind_r, 1) = P_tree(j, 2);
         
         %Find downstream pressure of child
         P_tree(ind_l, 2) = P_tree(j, 2) - dP_tree(ind_l);
         P_tree(ind_r, 2) = P_tree(j, 2) - dP_tree(ind_r);
         
         %Get thickness
         h_l = thickness(seg_size(ind_l, 3));
         h_r = thickness(seg_size(ind_r, 3));
         
         Sigma_tree(ind_l) = (P_tree(ind_l, 1) + P_tree(ind_l, 2)) / 2 * seg_size(ind_l,1) / h_l;
         Sigma_tree(ind_r) = (P_tree(ind_r, 1) + P_tree(ind_r, 2)) / 2 * seg_size(ind_r,1) / h_r;
     end
     if (seg_connectivity(j,1) == 0)
         %Find root stress
         Sigma_tree(j) = (P_tree(j, 1) + P_tree(j, 2)) / 2 * seg_size(j,1) / thickness(seg_size(j, 3));
     end
     
 end
 
 [seg_size_sort_ord, I] = sort(seg_size(:,3));
 %seg_size_sort = seg_size(I,:);
 
 for ord = 1:n_orders
     orig_ind_for_ord = I( seg_size_sort_ord == ord );
     n_seg_ord(ord) = length(orig_ind_for_ord);
     
     P_order = (P_tree(orig_ind_for_ord, 1) + P_tree(orig_ind_for_ord, 2)) / 2;
     P_seg_ord(ord, 1) = mean(P_order);
     P_seg_ord(ord, 2) = std(P_order);
     
     Q_order = Q_tree(orig_ind_for_ord);
     Q_seg_ord(ord, 1) = mean(Q_order);
     Q_seg_ord(ord, 2) = std(Q_order);
     
     Sigma_order = Sigma_tree(orig_ind_for_ord);
     Sigma_seg_ord(ord, 1) = mean(Sigma_order);
     Sigma_seg_ord(ord, 2) = std(Sigma_order);
     
     WSS_order = WSS_tree(orig_ind_for_ord);
     WSS_seg_ord(ord, 1) = mean(WSS_order);
     WSS_seg_ord(ord, 2) = std(WSS_order);
     
 end
 
 %Sorted by order
 hemo.P_order = P_seg_ord;
 hemo.Q_order = Q_seg_ord;
 hemo.Sigma_order = Sigma_seg_ord;
 hemo.WSS_order = WSS_seg_ord;
 
 %For the whole tree
 hemo.P_tree = (P_tree(:,1) + P_tree(:,2)) / 2;
 hemo.Q_tree = Q_tree;
 hemo.Sigma_tree = Sigma_tree;
 hemo.WSS_tree = WSS_tree;
 hemo.SegRad = seg_size(:,1);
 
end


function rel_vis=rel_viscosity(d)
%Pries et al Resistance to blood flow in microvessel in vivo.
%Assume Hd=0.45
rel_vis= (1+(6*exp(-0.0858*d)+3.2-2.44*exp(-0.06*d^0.645)-1)*(d/(d-1.1))^2)*(d/(d-1.1))^2;
end