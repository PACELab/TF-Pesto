% MILP for determining the placement and scheduling of ops.
% output: 
%a .mat of scheduling output; used for simulator;
%a .csv of scheduling output; for Ubaid experiment
%using v5solver
 
addpath ('D:\Program Files\IBM\ILOG\CPLEX_Studio1210\cplex\matlab\x64_win64')
strType = {'nmt-2lyr-1024-PCIE'};
[~,num_models] = size(strType);
pci = 1;
 
timelmt = 200;

%[# size of dag; size of largest node; longest length of merged supernode]
para_config = [700; 12000; 80];
[~,size_config] = size(para_config);

for fileseq = 1:1
    for cur_size = 1:size_config
        
fprintf('Model %s starts at %s\n', strType{fileseq},datestr(now,'HH:MM:SS.FFF'))
clearvars -except fileseq strType num_models  para_config size_config cur_size pci timelmt
%allowed_time = 12000; 
num_c2g_links = 2;
num_gpus = 2;

if pci
    load('speed_mat_pinned.mat');
    strPath = 'D:\GPU_Scheduling\logs\pcie';
    my_directory = 'D:\GPU_Scheduling\results\pci\coarsening\';
else    
    load('comm_nvlink02_pinned.mat');
    strPath = 'D:\GPU_Scheduling\logs\nvlink_02';
    my_directory = 'D:\GPU_Scheduling\results\NVLink_02\Coarsening\';
end

strName = {'edge-newk_cpu_ops_ids.json','edge-newk_gpu_ops_ids.json','edge-updated_execk_graph-final_ids.json','edge-newk_gpukernel_ops_ids.json','colocs.json','memory.json'};
 

% upper bound for total ops
total_expected_ops = 120000; 

% mapping from x1,x2, to ops_id position
map_ops = zeros(total_expected_ops,1);

coalloc = 1;
disp('Reading json log...')
tic;
cpu_struct = jsondecode(fileread(fullfile(strPath,strType{fileseq},strName{1})));
gpu_struct = jsondecode(fileread(fullfile(strPath,strType{fileseq},strName{2})));
graph_struct = jsondecode(fileread(fullfile(strPath,strType{fileseq},strName{3})));
kernel_struct = jsondecode(fileread(fullfile(strPath,strType{fileseq},strName{4})));

temp_filepath = fullfile(strPath,strType{fileseq},strName{5});
if isfile(temp_filepath)
coalloc_struct = jsondecode(fileread(temp_filepath));
else
    coalloc = 0;
end
memory_struct = jsondecode(fileread(fullfile(strPath,strType{fileseq},strName{6})));

% get ops id
names_cpu = fieldnames(cpu_struct);
cpu_total_ops = numel(names_cpu);
names_kernel = fieldnames(kernel_struct);
kernel_total_ops = numel(names_kernel);
names_gpu = fieldnames(gpu_struct);
gpu_total_ops = numel(names_gpu);



% comp = [type,median_processing_time]
% type: 1:cpu 2:gpu 3:kernel with a GPU exec 4. kernel without a GPU exec
size_node = cpu_total_ops+ gpu_total_ops+kernel_total_ops;
comp = zeros(size_node,2);

origin_order = zeros(size_node,1);
%disp('Fetching vertex information now...')

for i = 1:cpu_total_ops
    temp_idd = names_cpu{i}(2:end);  % remove the front "x"
    temp_loc  = str2double(temp_idd);
    origin_order(i) = temp_loc;
    comp(temp_loc,1) = 1;  %cpu
    comp(temp_loc,2) = cpu_struct.(names_cpu{i}).median;
end
for i = 1:kernel_total_ops
    temp_idd = names_kernel{i}(2:end);  % remove the front "x"
    temp_loc  = str2double(temp_idd);
    origin_order(cpu_total_ops+i) = temp_loc;
    comp(temp_loc,1) = 4;  % will be modified later if it has a gpu exec
    comp(temp_loc,2) = kernel_struct.(names_kernel{i}).median;
end
for i = 1:gpu_total_ops
    temp_idd = names_gpu{i}(2:end);
    temp_loc  = str2double(temp_idd);
    origin_order(cpu_total_ops+kernel_total_ops+i) = temp_loc;
    comp(temp_loc,1) = 2;  %gpu
    comp(temp_loc,2) = gpu_struct.(names_gpu{i}).median;
end

elapsed = toc;
fprintf('Vertex infomation: %g seconds \n', elapsed);

%disp('Fetching edge information now...')


%attached after current node (nodes after target node in a supernode)
attach_list = cell(60000,1);

names_graph = fieldnames(graph_struct);

% adjacency-lists [node array, size array]
prec_list = cell(size_node,2);
succ_list = cell(size_node,2);  

super_n = 1:size_node;  %point to its parent

% specific paring for kernel and its gpu exec
paring = zeros(size_node,1);

%logical for adj_topo
%adj_topo = false(size_node);
%adj_size = zeros(size_node);


for i = 1:numel(names_graph)
    temp_idd = names_graph{i}(2:end); %get rid of "x" again
    temp_loc  = str2double(temp_idd);

    current_node = names_graph{i};
    prece_nodes = fieldnames(graph_struct.(current_node));
    if ~isempty(prece_nodes)
        indeg = numel(prece_nodes);  %indegree
        for prece_entry = 1: indeg
            prece_id = prece_nodes{prece_entry}(2:end); 
            prece_loc  = str2double(prece_id);
            
            temp = fieldnames(graph_struct.(current_node).(prece_nodes{prece_entry}));
            temp_val = graph_struct.(current_node).(prece_nodes{prece_entry}).(temp{1}); % field includes dur,st and size
            
            if comp(temp_loc,1) == 2 % gpu ops, find j's precedence
                paring(prece_loc) = temp_loc;  % so 
                comp(prece_loc,1) = 3;  % kernel with a gpu exec
                comp(prece_loc,2) = 1;  % due to the overlap, changing the compute time to 1 for kernel with gpu
            end
            
            succ_list{prece_loc,1} = [succ_list{prece_loc,1} temp_loc];
            prec_list{temp_loc,1} = [prec_list{temp_loc,1} prece_loc];
          %  adj_topo(prece_loc,temp_loc) = 1;
            if isequal(temp_val,0)
                succ_list{prece_loc,2} = [succ_list{prece_loc,2} 0];
                prec_list{temp_loc,2} = [prec_list{temp_loc,2} 0];
            else
         %       adj_size(prece_loc,temp_loc) = temp_val.size(1); 
                succ_list{prece_loc,2} = [succ_list{prece_loc,2}  temp_val.size(1)];
                prec_list{temp_loc,2} = [prec_list{temp_loc,2} temp_val.size(1)];
            end                
        end
    end
end
origin_comp = comp;
elapsed = toc;
fprintf('Edge infomation: %g seconds \n', elapsed);

old_prec_list = prec_list;
old_succ_list = succ_list;

%recording colocs
if  coalloc
names_coalloc = fieldnames(coalloc_struct);
coalloc_total = numel(names_coalloc);
equal_colloc = cell(coalloc_total,1);

    for i = 1: coalloc_total
        t_length = numel(coalloc_struct.(names_coalloc{i}));
        temp_array = coalloc_struct.(names_coalloc{i});
        for j = 1: t_length
           tmp_str = str2double(temp_array{j});  
           if comp(tmp_str,1) > 1  %ignore cpu
           equal_colloc{i} = [equal_colloc{i} tmp_str];
           end
        end
    end

elapsed = toc;
fprintf('Parsing colocs: %g seconds \n', elapsed);
end

%step 1: merging kernel and its gpu exec
for i = 1:size_node
    gpu_id = paring(i);
    if gpu_id > 0  % kernel with a gpu exec        
        [super_n,comp,prec_list,succ_list,attach_list] = ops_merge_type(super_n,comp,prec_list,succ_list,attach_list,i,gpu_id);  

    end
end

elapsed = toc;
fprintf('Kernel GPU merging: %g seconds \n', elapsed);
fprintf('After kernel Coarsening: %u ops \n', nnz(comp(:,1)));

% %step 2: merging colocs as many as possible
% if  coalloc
%     for i = 1: coalloc_total        
%         target_array = unique(super_n(equal_colloc{i}));
%         level_ops = topo_level_v2(comp,prec_list,size_node);
%         tar_length  = length(target_array);
%         
%         tmp_matrix = sortrows([target_array' level_ops(target_array)],2,'descend');
%         for j = 1:tar_length-1
%             [super_n,comp,prec_list,succ_list,attach_list] = ops_merge_type(super_n,comp,prec_list,succ_list,attach_list,tmp_matrix(j+1,1),tmp_matrix(j,1));
%                      
%         end
%         
%         
%         
%     end
% 
% fprintf('After colocs Coarsening: %u ops \n', nnz(comp(:,1)));
% end

%step 3:regular merging
while nnz(comp(:,1))>para_config(1,cur_size)
good_val = nnz(comp(:,1));
level_ops = topo_level_v2(comp,prec_list,size_node);
marked = zeros(size_node,1);

edge_store =zeros(500000,3);
edge_ct = 1;
%generating all edges
for i = 1:size_node
   if comp(i,1)>0 
       suc_length =  length(succ_list{i,1});
       if suc_length>0
          edge_store(edge_ct:edge_ct+suc_length-1,1) = i;
          edge_store(edge_ct:edge_ct+suc_length-1,2) = succ_list{i,1}';
          edge_store(edge_ct:edge_ct+suc_length-1,3) = succ_list{i,2}';
          edge_ct = edge_ct + suc_length;
       end
   end
end
edge_ct = edge_ct-1;
edge_store = edge_store(1:edge_ct,:);

%sort based on size
[~,I] = sortrows(edge_store,3,'descend');

total_batch = zeros(5000,2);
row_ct = 1;
%merge if level_ops(successor) = level_ops(predecessor)+1;
for i = 1:edge_ct
   %start from largest comm
   cur_row = I(i);
   father = edge_store(cur_row,1);
   son = edge_store(cur_row,2);
   
   if marked(father)==1 || marked(son)==1 || (comp(father,1)==1 && comp(son,1)> 1) || (comp(son,1)==1 && comp(father,1)> 1) 
       continue
   end
   ops_sort = sort(level_ops(succ_list{father,1}));
       if  (length(prec_list{son,1})==1 || length(succ_list{father,1})==1 || level_ops(son) == level_ops(father)+1 || level_ops(son)<ops_sort(2) )  && comp(father,2)+ comp(son,2) < para_config(2,cur_size) && length(attach_list{father})<para_config(3,cur_size)
            marked(father) =1;
            marked(son) = 1;
            total_batch(row_ct,:) = [father son];
            row_ct = row_ct +1;
           pre_set = prec_list{son,1};
           pre_size = length(pre_set);
           for j = 1:pre_size
               parent = pre_set(j);
               if level_ops(parent) == level_ops(father)
                   marked(parent) = 1;
               end
           end
%             if row_ct>max(100,good_val*0.3)
%                 break
%             end
                
       end
end
if total_batch(1,1) == 0
    para_config(2,cur_size) = para_config(2,cur_size)+1000;
    para_config(3,cur_size) = para_config(3,cur_size)+1;
    continue
else
    [r,~] = size(total_batch);
    for z = 1:row_ct-1
%         a = total_batch(z,1);
%         b = total_batch(z,2);
%         if ~ismember(a,prec_list{b,1}) || ~ismember(b,succ_list{a,1})
%             disp('bugbug')
%             a 
%             b
%         end
        
        [super_n,comp,prec_list,succ_list,attach_list] = ops_merge_type(super_n,comp,prec_list,succ_list,attach_list,total_batch(z,1),total_batch(z,2));   
    end
end
end
para_config
fprintf('Coarsening graph has %u ops \n', nnz(comp(:,1)));


% extending the graph




% data to create
% a adjancy matrix (0-1matrix) for all ops
% a cpt matrix for all variables in the ILP
% a prec and succ list
xtotal = 20000;
cpt = zeros(10^7,5);
prec = cell(xtotal,1);
% succ = cell(xtotal,1);
adj = zeros(xtotal);

sep = zeros(7,1);



left_ops = nnz(comp(:,1));
l_map = zeros(size_node,1);

cnt = 1;
for i = 1:size_node
    if comp(i,1) > 0
       l_map(i) = cnt;
       cpt(cnt,1:2) = comp(i,:);
       cnt = cnt+1;
    end
end

% type: 1cpu 2gpu 4kernel 5c2g 6g2g
cnt = left_ops + 1;
for i = 1:size_node
   if ~isempty(succ_list{i,1})
       cur = l_map(i);
      cur_list = succ_list{i,1};
      q = length(cur_list);
      for j = 1:q
         child =  l_map(cur_list(j));
         if cpt(cur,1) == 1
             if cpt(child,1) == 1
             % cpu to cpu
             prec{child} = [prec{child} cur];
             adj(cur,child) = cur;
             else
                 % cpu to gpu, add one vertex
             cpt(cnt,1) = 5;
             cpt(cnt,2) = round(speed_mat(1,2)*succ_list{i,2}(j) + speed_mat(1,1));
             prec{child} = [prec{child} cnt];
             prec{cnt} = cur;
             adj(cur,cnt) = 1;
             adj(cnt,child) = 1;
             cnt = cnt + 1;
             end
         else
             if cpt(child,1) == 1
             % g2c
             cpt(cnt,1) = 5;
             cpt(cnt,2) = round(speed_mat(1,2)*succ_list{i,2}(j) + speed_mat(1,1));
             prec{child} = [prec{child} cnt];
             prec{cnt} = cur;
             adj(cur,cnt) = 1;
             adj(cnt,child) = 1;
             cnt = cnt + 1; 
             else
             % g2g
             cpt(cnt,1) = 6;
             cpt(cnt,2) = round(max(15,speed_mat(2,2)*succ_list{i,2}(j) + speed_mat(2,1)));
             cpt(cnt,3) = cur;
             cpt(cnt,4) = child;
             prec{child} = [prec{child} cnt];
             prec{cnt} = cur;
             adj(cur,cnt) = 1;
             adj(cnt,child) = 1;
             cnt = cnt + 1;     
             end    
         end
      end     
   end   
end
% total ops including added one
total_added_ops = cnt-1;
sep(1) = total_added_ops;
adj = adj(1:total_added_ops,1:total_added_ops);
prec = prec(1:total_added_ops);

% for i = 1:total_added_ops  %only for testing
%    if cpt(i,1)>=5
%        cpt(i,2)= cpt(i,2)*0.95;
%    end
%     
% end



% skip one for C_max
cnt = cnt + 1;
sep(2) = total_added_ops + 1;
%get p(x) as placement for kernel and gpu ops
for i = 1:total_added_ops
   if cpt(i,1) == 2 || cpt(i,1) == 4 % add a plcement variable
      cpt(cnt,1) = i;
      cpt(i,3) = cnt; %cross referencing
      cnt = cnt + 1;       
   end    
end
sep(3) = cnt -1;
% get z_ij for g2g communication
for i = (left_ops+1):total_added_ops
    if cpt(i,1) == 6
       cpt(cnt,1) = i;
       cpt(cnt,2) = cpt(i,3);
       cpt(cnt,3) = cpt(i,4);
       %(i,2) is processing time, 3 4 is prec and succ
       cpt(i,5) = cnt;
       cnt = cnt + 1;
    end
    
end
sep(4) = cnt -1;

%generate a starting point: x_0 by placing all jobs in gpu0
exp_G = digraph(adj);
topos = toposort(exp_G);

[exp_size,~] = size(adj);
x_0 = zeros(10^7,1);
cur_time = 0;
for i = 1:exp_size
    cur = topos(i);
    if cpt(cur,1) == 6
        x_0(cur) = x_0( cpt(cur,3));
    else
        x_0(cur) = cur_time + cpt(cur,2);
        cur_time = x_0(cur);
    end
end
x_0(exp_size+1) = max(x_0(1:exp_size));






% cpu congestion
[cntivity] = graphallshortestpaths(sparse(adj));
for i = 1:(total_added_ops-1)
   for j = (i+1): total_added_ops
       if cpt(i,1) == 1 && cpt(j,1) == 1 && cntivity(i,j) == inf && cntivity(j,i) == inf
       cpt(cnt,1:2) = [i j];
       if x_0(i) < x_0(j)
           x_0(cnt) = 1;
       end           
       cnt = cnt + 1;
       end
   end   
end
sep(5) = cnt -1;
%gpu congestion
for i = 1:(left_ops-1)
   for j = (i+1): left_ops
       if cpt(i,1) == 2 && cpt(j,1) == 2 && cntivity(i,j) == inf && cntivity(j,i) == inf
       cpt(cnt,1:2) = [i j];
       if x_0(i) < x_0(j)
           x_0(cnt) = 1;
       end         
       cnt = cnt + 1;
       end
   end   
end
sep(6) = cnt -1;
% communication congestion
for i = 1:(total_added_ops-1)
   for j = (i+1): total_added_ops
       if cpt(i,1) == 6 && cpt(j,1) == 6 && cntivity(i,j) == inf && cntivity(j,i) == inf
       cpt(cnt,1:2) = [i j];
       if x_0(i) < x_0(j)
           x_0(cnt) = 1;
       end    
       cnt = cnt + 1;
       end
   end   
end

cnt = cnt -1;
sep(7) = cnt;
cpt = cpt(1:cnt,1:end);
x_0 = x_0(1:cnt);






if coalloc
    for i = 1:coalloc_total
        equal_colloc{i} = unique(super_n(equal_colloc{i}));
        equal_colloc{i} = l_map(equal_colloc{i});
    end


[qqq,fval,exitflag,output]=v5_solver(prec,cpt,sep,x_0,timelmt,equal_colloc);
else
    [qqq,fval,exitflag,output]=v3_solver(prec,cpt,sep,x_0,timelmt);
end
qqq(sep(2)+1:end) = round(qqq(sep(2)+1:end));
%load('qqq.mat');

resu = qqq(1:left_ops);

[~,sorted]= sort(resu);

placem = zeros(size_node,2);
ctt = 1;
for i = 1:left_ops
    k = find(l_map == sorted(i),1);
    attach_list{k} = [k attach_list{k}];
    if origin_comp(k,1) == 1
        target_dev = 0;
    else
    target_dev = qqq(cpt(sorted(i),3));
    end
    for j = 1:length(attach_list{k})
            
            if origin_comp(attach_list{k}(j),1)== 1
                placem(attach_list{k}(j),1) = 1;

            elseif origin_comp(attach_list{k}(j),1)== 2
                placem(attach_list{k}(j),1) = 6+2*target_dev;
            else
                placem(attach_list{k}(j),1) = 2+2*target_dev;
            end                   
                placem(attach_list{k}(j),2) = ctt;
                ctt = ctt+1;
    end
        
end   
    

ops_id = cell(size_node,1);
final_output = zeros(size_node,2);
for i = 1:size_node
    asd = origin_order(i);
    final_output(i,:) = placem(asd,:);
    ops_id{i} = strcat('x',num2str(asd));
    
end




T = table(ops_id,final_output);
% my_directory = 'C:\Users\gabri\Documents\MATLAB\GPU_Scheduling\Critical_Path\NVlink_result\';  
temp_str = strcat(my_directory, num2str(fval),'-',strType{fileseq},'.csv');
writetable(T,temp_str,'Delimiter',',','WriteVariableNames',false);

temp_str2 = strcat(my_directory, num2str(fval),'-',strType{fileseq},'.mat');
save(temp_str2,'final_output');

    end
end
