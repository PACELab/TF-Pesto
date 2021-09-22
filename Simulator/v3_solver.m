% ILP information required

%1  standard ops (k,gpu,+added_nodes): completion time; (a pointer to the p(x)
%2  C_max  1;
%3  placement p(x) for kernel and gpu ops;
%4  z_ij() for g2g communication (prec,succ)
%5  delta_ij for cpu pair (a,b)
%6  theta_ij for gpu pair (a,b)
%7  omega_ij for g2g comm pair  (pair 1 pair 2)

%cpu ops: [1, pt];
%gpu ops: [2, pt, place]
%c2g comm: [5, pt];
%g2g comm [6 pt prec succ z(existence)]
function [x,fval,exitflag,output] = v3_solver(prec,cpt,sep,x_0,timelmt)


M_s = 2*10^6;
M_l = 4*10^6;

total_var = sep(7);
total_ops = sep(1);

%1---precedence constraints
B_neq1 = zeros(150000,1);
total_var_row1 = zeros(600000,3);
row_neq1 = 1;
row_sparse_idx = 1;
for j = 1: total_ops
        if ~isempty(prec{j}) 
            % c_i + p_j - c_j <= 0
            for i = 1:numel(prec{j}) %i -> j
                %[row column val]
                total_var_row1(row_sparse_idx,:) = [row_neq1 prec{j}(i) 1];
                row_sparse_idx = row_sparse_idx+1;
                
                    if cpt(j,1) == 6   % g2g nodes
                        % p_jz_j
                        total_var_row1(row_sparse_idx,:) = [row_neq1 cpt(j,5) cpt(j,2)];
                        row_sparse_idx = row_sparse_idx+1;
                    else
                        B_neq1(row_neq1) = -cpt(j,2);  %-p_j                      
                    end
                    % -c_j
                total_var_row1(row_sparse_idx,:) = [row_neq1 j -1];    
                row_sparse_idx = row_sparse_idx+1;     

                row_neq1 = row_neq1+1;
            end
        end
end
total_var_row1(row_sparse_idx:end,:) = []; 
B_neq1(row_neq1:end) = [];
SP_neq1 = sparse(total_var_row1(:,1),total_var_row1(:,2),total_var_row1(:,3),row_neq1-1,total_var);



%2---makespan
% c_j - C <= 0
total_var_row2 = zeros(2*total_ops,3);
for i = 1:total_ops
    total_var_row2(2*i-1,:) = [i i 1];
    total_var_row2(2*i,:) = [i sep(2) -1];    
end
SP_neq2 = sparse(total_var_row2(:,1),total_var_row2(:,2),total_var_row2(:,3),total_ops,total_var);
B_neq2 = zeros(1,total_ops);

%3---no earlier than completion time
% -c_i <= -p_i
total_var_row3 = zeros(total_ops,3);
B_neq3 = zeros(total_ops,1);
for i = 1: total_ops
    total_var_row3(i,:) = [i i -1];
    B_neq3(i) = -cpt(i,2);    
end
SP_neq3 = sparse(total_var_row3(:,1),total_var_row3(:,2),total_var_row3(:,3),total_ops,total_var);

%4---xor condition for existence of comm
% z_ij =1 iff y_i and y_j differ

total_var_row4 = zeros(600000,3);
row_neq4 = 1;
temp_row_idx = 1;
% z: cpt(i,5), x:cpt(cpt(i,3),3) y: cpt(cpt(i,4),3))
for i = 1:total_ops
    if cpt(i,1) == 6
        % z- x-y <=0
        total_var_row4(temp_row_idx,:) = [row_neq4 cpt(i,5) 1];
        total_var_row4(temp_row_idx+1,:) = [row_neq4 cpt(cpt(i,3),3) -1];
        total_var_row4(temp_row_idx+2,:) = [row_neq4 cpt(cpt(i,4),3) -1];
        row_neq4 = row_neq4+ 1;
        temp_row_idx = temp_row_idx+3;
        % -z+x-y<=0
        total_var_row4(temp_row_idx,:) = [row_neq4 cpt(i,5) -1];
        total_var_row4(temp_row_idx+1,:) = [row_neq4 cpt(cpt(i,3),3) 1];
        total_var_row4(temp_row_idx+2,:) = [row_neq4 cpt(cpt(i,4),3) -1];
        row_neq4 = row_neq4+ 1;
        temp_row_idx = temp_row_idx+3;        
        % -z-x+y <=0
        total_var_row4(temp_row_idx,:) = [row_neq4 cpt(i,5) -1];
        total_var_row4(temp_row_idx+1,:) = [row_neq4 cpt(cpt(i,3),3) -1];
        total_var_row4(temp_row_idx+2,:) = [row_neq4 cpt(cpt(i,4),3) 1];
        row_neq4 = row_neq4+ 1;
        temp_row_idx = temp_row_idx+3;  
        % z+x+y <=2
        total_var_row4(temp_row_idx,:) = [row_neq4 cpt(i,5) 1];
        total_var_row4(temp_row_idx+1,:) = [row_neq4 cpt(cpt(i,3),3) 1];
        total_var_row4(temp_row_idx+2,:) = [row_neq4 cpt(cpt(i,4),3) 1];
        row_neq4 = row_neq4+ 1;
        temp_row_idx = temp_row_idx+3; 
    end 
end
total_var_row4(temp_row_idx:end,:) = []; 
B_neq4 = repmat([0 0 0 2],1,(row_neq4-1)/4);
SP_neq4 = sparse(total_var_row4(:,1),total_var_row4(:,2),total_var_row4(:,3),row_neq4-1,total_var);




%5--- congestion on cpu
st = sep(4)+1;
ed = sep(5);
total_var_row5 = zeros(600000,3);
B_neq5 = zeros(2*(ed-st+1),1);
row_neq5 = 1;
temp_row_idx = 1;
for i = st:ed
   jobA = cpt(i,1);
   jobB = cpt(i,2);
   % -c_i + c_j - M*delta_ij <= -p_i
   total_var_row5(temp_row_idx,:) = [row_neq5 jobA -1];
   total_var_row5(temp_row_idx+1,:) = [row_neq5 jobB 1];
   total_var_row5(temp_row_idx+2,:) = [row_neq5 i -M_s];
   B_neq5(row_neq5) = -cpt(jobA,2);
   row_neq5 = row_neq5 + 1;
   temp_row_idx = temp_row_idx + 3;
   % c_i -c_j + M*delta_ij <= M-P_j
   total_var_row5(temp_row_idx,:) = [row_neq5 jobA 1];
   total_var_row5(temp_row_idx+1,:) = [row_neq5 jobB -1];
   total_var_row5(temp_row_idx+2,:) = [row_neq5 i M_s];
   B_neq5(row_neq5) = -cpt(jobB,2)+ M_s;
   row_neq5 = row_neq5+ 1;
   temp_row_idx = temp_row_idx + 3;  
end
total_var_row5(temp_row_idx:end,:) = []; 
SP_neq5 = sparse(total_var_row5(:,1),total_var_row5(:,2),total_var_row5(:,3),row_neq5-1,total_var);



%%6 ---congestion on gpu
st = sep(5)+1;
ed = sep(6);
total_var_row6 = zeros(30000000,3);
B_neq6 = zeros(4*(ed-st+1),1);
row_neq6 = 1;
temp_row_idx = 1;
for i = st:ed
   jobi = cpt(i,1);
   jobj = cpt(i,2);
    % -c_i + c_j - M*delta_ij + M_L*I(i) + M_L*I(j) <= -p_i + 2*M*L
   total_var_row6(temp_row_idx,:) = [row_neq6 jobi -1];
   total_var_row6(temp_row_idx+1,:) = [row_neq6 jobj 1];
   total_var_row6(temp_row_idx+2,:) = [row_neq6 i -M_s];
   total_var_row6(temp_row_idx+3,:) = [row_neq6 cpt(jobi,3) M_l];
   total_var_row6(temp_row_idx+4,:) = [row_neq6 cpt(jobj,3) M_l];
   B_neq6(row_neq6) = -cpt(jobi,2)+ 2*M_l;
   row_neq6 = row_neq6+1;
   temp_row_idx = temp_row_idx + 5;
   % c_i -c_j + M*delta_ij + M_l*I(i) + M_l*I(j) <= M-P_j + 2*M_l
   total_var_row6(temp_row_idx,:) = [row_neq6 jobi 1];
   total_var_row6(temp_row_idx+1,:) = [row_neq6 jobj -1];
   total_var_row6(temp_row_idx+2,:) = [row_neq6 i M_s];
   total_var_row6(temp_row_idx+3,:) = [row_neq6 cpt(jobi,3) M_l];
   total_var_row6(temp_row_idx+4,:) = [row_neq6 cpt(jobj,3) M_l];
   B_neq6(row_neq6) = M_s-cpt(jobj,2)+ 2*M_l;
   row_neq6 = row_neq6+1;
   temp_row_idx = temp_row_idx + 5;   
   % -c_i + c_j - M*delta_ij-M_l*I(i) - M_l*I(j) <= -p_i
   total_var_row6(temp_row_idx,:) = [row_neq6 jobi -1];
   total_var_row6(temp_row_idx+1,:) = [row_neq6 jobj 1];
   total_var_row6(temp_row_idx+2,:) = [row_neq6 i -M_s];
   total_var_row6(temp_row_idx+3,:) = [row_neq6 cpt(jobi,3) -M_l];
   total_var_row6(temp_row_idx+4,:) = [row_neq6 cpt(jobj,3) -M_l];
   B_neq6(row_neq6) = -cpt(jobi,2);
   row_neq6 = row_neq6+1;
   temp_row_idx = temp_row_idx + 5;   
   % c_i -c_j + M*delta_ij -M_l*I(i) - M_l*I(j) <= M-P_j
   total_var_row6(temp_row_idx,:) = [row_neq6 jobi 1];
   total_var_row6(temp_row_idx+1,:) = [row_neq6 jobj -1];
   total_var_row6(temp_row_idx+2,:) = [row_neq6 i M_s];
   total_var_row6(temp_row_idx+3,:) = [row_neq6 cpt(jobi,3) -M_l];
   total_var_row6(temp_row_idx+4,:) = [row_neq6 cpt(jobj,3) -M_l];
   B_neq6(row_neq6) = M_s-cpt(jobj,2);
   row_neq6 = row_neq6+1;
   temp_row_idx = temp_row_idx + 5;       
end
total_var_row6(temp_row_idx:end,:) = []; 
SP_neq6 = sparse(total_var_row6(:,1),total_var_row6(:,2),total_var_row6(:,3),row_neq6-1,total_var);



% congestion on comm
%%6 ---congestion on gpu
st = sep(6)+1;
ed = sep(7);
total_var_row7 = zeros(300000000,3);
B_neq7 = zeros(2*(ed-st+1),1);
row_neq7 = 1;
temp_row_idx = 1;
for i = st:ed
   commi = cpt(i,1);
   commj = cpt(i,2);
   
   joba = cpt(cpt(commi,3),3); %prec's placement index
   jobb = cpt(cpt(commi,4),3);
   jobc = cpt(cpt(commj,3),3);
   jobd = cpt(cpt(commj,4),3);
   % (c_i-p_i)+Mdelta_ij -c_j - M(a+c-b-d-2)>=0
   total_var_row7(temp_row_idx,:) = [row_neq7 commi -1];
   total_var_row7(temp_row_idx+1,:) = [row_neq7 commj 1];
   total_var_row7(temp_row_idx+2,:) = [row_neq7 i -M_s];
   total_var_row7(temp_row_idx+3,:) = [row_neq7 joba M_l];
   total_var_row7(temp_row_idx+4,:) = [row_neq7 jobb -M_l];
   total_var_row7(temp_row_idx+5,:) = [row_neq7 jobc M_l];
   total_var_row7(temp_row_idx+6,:) = [row_neq7 jobd -M_l];
   B_neq7(row_neq7) = -cpt(commi,2)+ 2*M_l;
   row_neq7 = row_neq7+1;
   temp_row_idx = temp_row_idx + 7;  
    
   total_var_row7(temp_row_idx,:) = [row_neq7 commi 1];
   total_var_row7(temp_row_idx+1,:) = [row_neq7 commj -1];
   total_var_row7(temp_row_idx+2,:) = [row_neq7 i M_s];
   total_var_row7(temp_row_idx+3,:) = [row_neq7 joba M_l];
   total_var_row7(temp_row_idx+4,:) = [row_neq7 jobb -M_l];
   total_var_row7(temp_row_idx+5,:) = [row_neq7 jobc M_l];
   total_var_row7(temp_row_idx+6,:) = [row_neq7 jobd -M_l];
   B_neq7(row_neq7) = M_s-cpt(commj,2)+ 2*M_l;
   row_neq7 = row_neq7+1;
   temp_row_idx = temp_row_idx + 7;  
   
   %(c_i-p_i)+Mdelta_ij -c_j - M(b+d-a-c-2)>=0
   total_var_row7(temp_row_idx,:) = [row_neq7 commi -1];
   total_var_row7(temp_row_idx+1,:) = [row_neq7 commj 1];
   total_var_row7(temp_row_idx+2,:) = [row_neq7 i -M_s];
   total_var_row7(temp_row_idx+3,:) = [row_neq7 joba -M_l];
   total_var_row7(temp_row_idx+4,:) = [row_neq7 jobb M_l];
   total_var_row7(temp_row_idx+5,:) = [row_neq7 jobc -M_l];
   total_var_row7(temp_row_idx+6,:) = [row_neq7 jobd M_l];
   B_neq7(row_neq7) = -cpt(commi,2)+ 2*M_l;
   row_neq7 = row_neq7+1;
   temp_row_idx = temp_row_idx + 7;  
   
   total_var_row7(temp_row_idx,:) = [row_neq7 commi 1];
   total_var_row7(temp_row_idx+1,:) = [row_neq7 commj -1];
   total_var_row7(temp_row_idx+2,:) = [row_neq7 i M_s];
   total_var_row7(temp_row_idx+3,:) = [row_neq7 joba -M_l];
   total_var_row7(temp_row_idx+4,:) = [row_neq7 jobb M_l];
   total_var_row7(temp_row_idx+5,:) = [row_neq7 jobc -M_l];
   total_var_row7(temp_row_idx+6,:) = [row_neq7 jobd M_l];
   B_neq7(row_neq7) = M_s-cpt(commj,2)+ 2*M_l;
   row_neq7 = row_neq7+1;
   temp_row_idx = temp_row_idx + 7;    
   
end
total_var_row7(temp_row_idx:end,:) = []; 
SP_neq7 = sparse(total_var_row7(:,1),total_var_row7(:,2),total_var_row7(:,3),row_neq7-1,total_var);

% extra

ctype = [];
for ttt = 1: sep(2)
    ctype = [ctype; 'C'];
end
for ttt = (sep(2)+1):sep(7)
    ctype = [ctype; 'B'];
end
ctype = ctype';


f = zeros(1,sep(7)); 
f(sep(2)) = 1;

lb = zeros(1,total_var);
ub = [Inf*ones(1,sep(2)) ones(1,sep(7)-sep(2))];

 SP_left_set1 = [ SP_neq1; SP_neq2; SP_neq3; SP_neq4; SP_neq5; SP_neq6;  SP_neq7];
 B_right_set1 = [ B_neq1' B_neq2 B_neq3' B_neq4 B_neq5' B_neq6' B_neq7'];

% SP_left_set1 = [SP_neq1; SP_neq2; SP_neq3; SP_neq4; SP_neq5; SP_neq6; SP_neq7 ];
% B_right_set1 = [B_neq1' B_neq2 B_neq3' B_neq4 B_neq5' B_neq6' B_neq7'];

opt = cplexoptimset('cplex');
opt.mip.tolerances.mipgap = 0.002;
opt.mip.strategy.probe = 3;
%opt.mip.display = 3;
 opt.Display = 'iter';
%opt.Algorithm = 'dual';
%opt.emphasis.mip = 1;

opt.timelimit = timelmt;
 [x,fval,exitflag,output] = cplexmilp(f,SP_left_set1,B_right_set1',[],[],[],[],[],lb,ub,ctype,x_0,opt);

end




