

function resistance=find_opt_morphometric_tree(order,order_remodeling,alpha,beta,gen)

% Create matrices to store the morphometric data from Huang et al. in CGS units.
global connectivity length diameter viscosity remainder total_elements_created

%Initialize vessel order metrics
con_dim = size(connectivity);
n_orders = con_dim(1); %Number of vessel orders, vessels decrease in size with decreasing order

% reduce number of child elements with order_remodeling to 1 to model pruning  
for i=order_remodeling:-1:1
   connectivity(i,:)=connectivity(i,:)*alpha;
end

for i=order_remodeling:-1:1
   diameter(i)=diameter(i)*beta;
end

remainder = zeros(n_orders, n_orders);

total_elements_created = 0;

plasma_viscosity=0.0124; % Poise

density=1.06; % g/cm^3

%generate morphometric tree for a given order
if (gen == 1)
    generate_morphometric_tree(order);
end
load morphometric_tree
%      clculate resistance
number_seg=size(seg_connectivity,1);
seg_input_resistance=zeros(number_seg,1);

for i=number_seg:-1:1
    %Update seg_size for the new diameter array
    if (gen == 0)
            seg_size(i,1) = diameter(seg_size(i,3)) / 2;
    end
    
    if seg_connectivity(i,2)==0 && seg_connectivity(i,3)==0
        d=seg_size(i,1)*2*10000;%convert to micron;
        viscosity=rel_viscosity(d)*plasma_viscosity;
        seg_input_resistance(i)=8*viscosity*seg_size(i,2)/(pi*(seg_size(i,1))^4);
    end
    if seg_connectivity(i,2)==-1 && seg_connectivity(i,3)~=-1
        %   disp('a seg has not been connected to its child element'); %one of child input resistance will be zero treat it like a serial circuit
        d=seg_size(i,1)*2*10000;%convert to micron;
        viscosity=rel_viscosity(d)*plasma_viscosity;
        seg_input_resistance(i)=8*viscosity*seg_size(i,2)/(pi*(seg_size(i,1))^4)+seg_input_resistance(seg_connectivity(i,3));
    end
    if seg_connectivity(i,3)==-1 && seg_connectivity(i,2)~=-1
        %  disp('a seg has not been connected to its child element'); %one of child input resistance will be zero treat it like a serial circuit
        d=seg_size(i,1)*2*10000;%convert to micron;
        viscosity=rel_viscosity(d)*plasma_viscosity;
        seg_input_resistance(i)=8*viscosity*seg_size(i,2)/(pi*(seg_size(i,1))^4)+seg_input_resistance(seg_connectivity(i,2));
    end
    if seg_connectivity(i,2)>0 && seg_connectivity(i,3)>0
        d=seg_size(i,1)*2*10000;%convert to micron;
        viscosity=rel_viscosity(d)*plasma_viscosity;
        seg_input_resistance(i)=8*viscosity*seg_size(i,2)/(pi*(seg_size(i,1))^4)+1/(1/seg_input_resistance(seg_connectivity(i,2)) + 1/seg_input_resistance(seg_connectivity(i,3)));
    end
end

number_seg=size(seg_connectivity,1);
resistance=seg_input_resistance(1);


save morphometric_tree seg_connectivity parent_seg_elem seg_size seg_input_resistance
end

function rel_vis=rel_viscosity(d)
%Pries et al Resistance to blood flow in microvessel in vivo.
%Assume Hd=0.45
rel_vis=(1+(6*exp(-0.0858*d)+3.2-2.44*exp(-0.06*d^0.645)-1)*(d/(d-1.1))^2)*(d/(d-1.1))^2;
end
