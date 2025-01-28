% This function creates the children of a parent element of a given order using the data
% of Jiang et al.  
% All calculations are performed with CGS (centimeters, grams, and seconds) units.  
% Create tree directly without calling recursive functions and lable
% segment and elements number and create a connectivity matrix 

% Z. L. Jiang, G. S. Kassab, and Y. C. Fung, "Diameter-defined Strahler system and
% connectivity matrix of the pulmonary arterial tree," J Appl Physiol, v. 76, p. 882-892,
% 1994.

function generate_morphometric_trees4(parent_order)
if parent_order <= 15
  N=30000000;
end
% if parent_order >13 
% N=15000000;
% end
if parent_order >11 && parent_order<= 13
N=5000000;
end
if parent_order >7 && parent_order <=11
N=400000;
end
if parent_order <=7 
N=10000;
end
% List the global variables so that they can be accessed.
global connectivity length diameter viscosity remainder total_elements_created

% From the study of Poiseuille flow, we have the following formula for
% the resistance of each segment.
%       resistance = 8 * dynamic viscosity * length of the segment / ( pi * radius^4 )
% The resistances of the segments and downstream elements must be combined to obtain the
% resistance at the inlet node of the parent element.  The resistance of each child
% element must be found to compute the resistance of the parent element.  This is done
% by calling the calculate_resistance_HW1_soln function with the order of the child
% element as input.  We also use the following results from analagous resistive
% electrical circuits.  The total resistance R_t of two resistors with resistance values
% R_1 and R_2 in series is their sum, R_t=R_1+R_2.  The total resistance of two resistors
% in parallel is R_t=1/(1/R_1+1/R_2).

% Calculate the number of child elements of each order by rounding the sum of the mean
% value in table 4 and the current remainder to the nearest whole number.  Since the
% data on the number of pruned elements is omitted from the paper, their contribution to
% the morphometry-based tree is neglected.  Add the number of child elements to the
% running total of elements created.  Update the remainders based on the number of child
% elements created.

children=connectivity(:,parent_order)+remainder(:,parent_order);
% the use of remainder will result in a difference in final total
% resistance because the tree geometry is slightly different from that
% created by the original recursive function because the additional
% segments from the remainder are added at different locations i.e. the local geometry is denpendent on the order of creating elements   
remainder(:,parent_order)=children-round(children); 
children=round(children);
total_elements_created=total_elements_created+children;
number_of_children=sum(children);

% Create a list of the orders of the child elements from smallest to largest.
child_orders=[];
for order=1:parent_order
    child_orders=[child_orders; order*ones(children(order),1)];
    
end

current_working_elem=1;
parent_seg_elem=zeros(N,3);
parent_seg_elem(1,:)=[0,0,parent_order]; % [ seg number that connects the inlet of the current elem, element position (left 2 or right 3, no parent segment connected 0) and current elem order; nothing is connected to the fist element inlet;
max_elem=1;
% first segment connectivity, [parent_seg_index, left_child_seg_index, right_child_seg_index], -1 indicates children element has not been
% connected
%seg_connectivity=[];

%seg_size;
max_seg=0;
seg_connectivity=zeros(N,3);
seg_size=zeros(N,3);

if number_of_children <2
    
    % the root element (large order) should at least 2 children
    disp('children for root < 2, stop');
    return
end


if number_of_children>=2
    for i=1:number_of_children-1
        if i<number_of_children-1
            seg_connectivity(max_seg+1,:)=[i-1 -1 i+1]; % split element into number_of_children-1 pieces and assign segment number
            
        else
            seg_connectivity(max_seg+1,:)=[i-1 -1 -1] ;
            
        end
        seg_size(max_seg+1,:)=[diameter(parent_order)/2, length(parent_order)/(number_of_children-1), parent_order ];
        max_seg=max_seg+1;
    end
end

if child_orders(end)==parent_order&size(child_orders,1)~=2  %reorder child elements
    child_orders=[child_orders(1:end-3); child_orders(end);
        child_orders(end-2:end-1)];
end

% add last two child elems to the outlet of the current elem
parent_seg_elem(max_elem+1:max_elem+2,:)=[max_seg, 2, child_orders(end);max_seg,3, child_orders(end-1)];
max_elem=max_elem+2;

child_orders=child_orders(1:end-2);
% For each remaining child element, add to the interior segment's left side and record the parent segment for this child element and child elem order
for i=1:size(child_orders,1)
    parent_seg_elem(max_elem+1,:)=[max_seg-i,2, child_orders(end+1-i)];
    max_elem=max_elem+1;
end

current_working_elem=current_working_elem+1;
while current_working_elem <=max_elem
    
    par_index=parent_seg_elem(current_working_elem,1);
    side=parent_seg_elem(current_working_elem,2);
    current_order=parent_seg_elem(current_working_elem,3);
    %      first_seg_index=size(seg_connectivity,1)+1;
    first_seg_index=max_seg+1;
    children=connectivity(:,current_order)+remainder(:,current_order);
    remainder(:,current_order)=children-round(children);
    children=round(children);
    total_elements_created=total_elements_created+children;
    number_of_children=sum(children);
    
    % Create a list of the orders of the child elements from smallest to largest.
    child_orders=[];
    for order=1:parent_order
        child_orders=[child_orders; order*ones(children(order),1)];
    end
    
    
    if number_of_children==0
        seg_connectivity(par_index,side)=first_seg_index;
        seg_connectivity(max_seg+1,:)=[par_index 0 0];
        seg_size(max_seg+1,:)=[diameter(current_order)/2  length(current_order) current_order];
        max_seg=max_seg+1;
    end
    
    
    
    % There is one order 1 child of an order 1 parent.
    % split the element into 2 and add the child element to the mid point of parent element
    if number_of_children==1 && child_orders(1)==1 && current_order==1
        seg_connectivity(par_index,side)=first_seg_index;
        seg_connectivity(max_seg+1:max_seg+2,:)=[par_index -1 first_seg_index+1 ; first_seg_index 0 0];
        seg_size(max_seg+1:max_seg+2,:)=[diameter(1)/2 length(1)/2 1; diameter(1)/2 length(1)/2 1];
        max_seg=max_seg+2;
        parent_seg_elem(max_elem+1,:)=[first_seg_index, 2, 1];
        max_elem=max_elem+1;
    elseif number_of_children==1 && child_orders(1)~=1 && current_order~=1
        % add to the middle
        seg_connectivity(par_index,side)=first_seg_index;
        seg_connectivity(max_seg+1:max_seg+2,:)=[par_index -1 first_seg_index+1; first_seg_index 0 0];
        seg_size(max_seg+1:max_seg+2,:)=[diameter(current_order)/2 length(current_order)/2 current_order; diameter(current_order)/2 length(current_order)/2 current_order];
        max_seg=max_seg+2;
        parent_seg_elem(max_elem+1,:)=[first_seg_index, 2, child_orders(1) ];
        max_elem=max_elem+1;
        %      disp('number_of_children=1 child_order>1 current_order>1')
        %       child_orders(1)
        %       current_order
    elseif number_of_children==1
        % add to the middle
        seg_connectivity(par_index,side)=first_seg_index;
        seg_connectivity(max_seg+1:max_seg+2,:)=[par_index -1 first_seg_index+1; first_seg_index 0 0];
        seg_size(max_seg+1:max_seg+2,:)=[diameter(current_order)/2 length(current_order)/2 current_order; diameter(current_order)/2 length(current_order)/2 current_order];
        max_seg=max_seg+2;
        parent_seg_elem(max_elem+1,:)=[first_seg_index, 2, child_orders(1) ];
        max_elem=max_elem+1;
    end
    
    if number_of_children>=2
        
        seg_size(max_seg+1:max_seg+number_of_children-1,:)=[ones(number_of_children-1,1)*(diameter(current_order)/2), ones(number_of_children-1,1)*length(current_order)/(number_of_children-1), ones(number_of_children-1,1)*current_order];
        if child_orders(end)==parent_order&size(child_orders,1)~=2
            child_orders=[child_orders(1:end-3); child_orders(end);
                child_orders(end-2:end-1)];
        end
        
        % assign seg number for the first piece of the current element and
        % connect to the parent segment
        
        seg_connectivity(par_index,side)=first_seg_index;
        
        
        if number_of_children ==2
            seg_connectivity(max_seg+1,:)=[par_index -1 -1];
            max_seg=max_seg+1;
        else
            for i=1:number_of_children-1
                if i==1
                    seg_connectivity(max_seg+1,:)=[par_index -1 first_seg_index+1];
                end
                if i>1 && i<number_of_children-1
                    seg_connectivity(max_seg+1,:)=[first_seg_index+i-2 -1 first_seg_index+i]; % split element into number_of_children-1 pieces and assign segment number
                end
                if i==number_of_children-1
                    seg_connectivity(max_seg+1,:)=[first_seg_index+i-2 -1 -1];
                end
                max_seg=max_seg+1;
            end
        end
        
        
        % add last two child elems to the outlet of the current elem
        parent_seg_elem(max_elem+1:max_elem+2,:)=[max_seg 2 child_orders(end) ;max_seg 3 child_orders(end-1)];
        max_elem=max_elem+2;
        
        child_orders=child_orders(1:end-2);
        % For each remaining child element, add to the interior segment's left side and record the parent segment for this child element and child elem order
        
        
        for i=1:size(child_orders,1)
            parent_seg_elem(max_elem+1,:)=[max_seg-i, 2, child_orders(end+1-i)];
            max_elem=max_elem+1;
        end
        
        
    end
    current_working_elem=current_working_elem+1;
    
end

%save morphometric tree
seg_connectivity=seg_connectivity(1:max_seg,:);
parent_seg_elem=parent_seg_elem(1:max_elem,:);
seg_size=seg_size(1:max_seg,:);
save morphometric_tree seg_connectivity parent_seg_elem seg_size
   
   
   
   
