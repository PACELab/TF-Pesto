% Operation Merge

function [super_n,comp,prec_list,succ_list,attach_list] = ops_merge_type(super_n,comp,prec_list,succ_list,attach_list,parent,child)


% attach gpu to kernel
succ_pos = find(succ_list{parent,1} == child,1); %given small size, find should be fast
% trans_size = succ_list{parent,2}(succ_pos);
attach_list{parent} = [attach_list{parent} child attach_list{child}];
super_n(child) = parent;
super_n(attach_list{child}) = parent;
attach_list{child} = [];
% clean comp

comp(parent,2) = comp(parent,2) + comp(child,2);

%change type to child, cpu ops is already considered in the batch process
if comp(parent,1) ~= 2 
comp(parent,1) = comp(child,1);
end
comp(child,:) = [0 0];


%cut connection from parent to child
succ_list{parent,1}(succ_pos) = [];
succ_list{parent,2}(succ_pos) = [];

% add parent's grandchild to succ_list{parent}
if ~isempty(succ_list{child,1})
    sssc =  succ_list{child,1};
    w = length(sssc);
    for p = 1:w
        [Lia,Locb] =ismember(sssc(p),succ_list{parent,1});
       if Lia %already in
           succ_list{parent,2}(Locb) = succ_list{parent,2}(Locb) + succ_list{child,2}(p);
       else
           succ_list{parent,1} = [succ_list{parent,1} succ_list{child,1}(p)];
           succ_list{parent,2} = [succ_list{parent,2} succ_list{child,2}(p)];
       end
    end
    
% change the predecessor of all grandchildren    
    q = length(succ_list{child,1});
    for j = 1: q
        succ_ops = succ_list{child,1}(j);
        loc = prec_list{succ_ops,1} == child;
        [Lia,Locb] =ismember(parent,prec_list{succ_ops,1});
        if Lia  %already in
            prec_list{succ_ops,1}(loc) = [];  % remove child
            prec_list{succ_ops,2}(Locb) = prec_list{succ_ops,2}(Locb) + prec_list{succ_ops,2}(loc) ;
            prec_list{succ_ops,2}(loc) = [];
        else   % replace position where child is         
            prec_list{succ_ops,1}(loc) = parent;
            prec_list{succ_ops,2}(loc) = prec_list{succ_ops,2}(loc);
        end
    end
end



%also add all predecessor of child to the predecessor of parent
if length(prec_list{child,1}) > 1
    sssc =  prec_list{child,1};
    w = length(sssc);
    for p = 1:w
        if sssc(p) ~= parent
            [Lia,Locb] =ismember(sssc(p),prec_list{parent,1});
            if Lia  %already in
                prec_list{parent,2}(Locb) = prec_list{parent,2}(Locb)+ prec_list{child,2}(p);
                scp = find(succ_list{sssc(p),1} == parent);
                succ_list{sssc(p),2}(scp) = succ_list{sssc(p),2}(scp) + prec_list{child,2}(p);
            else
                prec_list{parent,1} = [prec_list{parent,1} sssc(p)];
                prec_list{parent,2} = [prec_list{parent,2} prec_list{child,2}(p)];
                succ_list{sssc(p),1} = [succ_list{sssc(p),1} parent];
                succ_list{sssc(p),2} = [succ_list{sssc(p),2} prec_list{child,2}(p)];
                              
            end
            
            remov = find(succ_list{sssc(p),1} == child);
            if remov > 0
            succ_list{sssc(p),1}(remov) = [];            
            succ_list{sssc(p),2}(remov) = [];
            end
        end
    end
           
    
end





prec_list{child,1} = [];
prec_list{child,2} = [];
succ_list{child,1} = [];
succ_list{child,2} = [];





end

