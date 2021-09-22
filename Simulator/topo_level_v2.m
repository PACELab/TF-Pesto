%ops in the same batch belongs to the same level
% root level = 1
function [level_ops] = topo_level_v2(comp,prec_list,size_node)
seq = [];
level_ops = zeros(size_node,1);
for i = 1:size_node
   if comp(i,1) > 0 && isempty(prec_list{i,1})   %head node
       seq = [seq i];
       level_ops(i) = 1;
   end 
   
end
marked = zeros(size_node,1);
marked(seq) = 1;

level_count = 2;
while ~isempty(seq)
    temp_seq = seq;
    seq = [];
    for i = 1:size_node
        if comp(i,1) > 0 && marked(i) == 0
        prec_list{i,1} = prec_list{i,1}(~ismember(prec_list{i,1},temp_seq));
        if isempty( prec_list{i,1}) 
            marked(i) = 1;
            level_ops(i) = level_count;
            seq = [seq i];
        end
        end
    end
    level_count = level_count +1;
end
end

