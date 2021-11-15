function full_M_all_lumped=ExtractMass()
load('fem_node_element.mat');load('alldata_Force_Displacement.mat');
shell3P = allshell(1:66,3:5);
shell4P = allshell(67:end,3:6);
newshell3P = zeros(66,3);
for i=1:66
    for j=1:3
        newshell3P(i,j) = find(node_num == shell3P(i,j));
    end
end
newshell4P = zeros(size(shell4P,1),4);
for i=1:size(shell4P,1)
    for j=1:4
        newshell4P(i,j) = find(node_num == shell4P(i,j));
    end
end
%calculate matrix M and F
M_all_1 = MassAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362);
M_all_2 = MassAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362);
M_all = M_all_1 + M_all_2;
%spy(M_all,'k')
full_M_all_original = kron(M_all,eye(3))*2.78*1.6*0.000001;
%spy(full_M_all_original,'k')
%lumped mass matrix
tem_matrix = sum(full_M_all_original);
full_M_all_lumped = diag(tem_matrix);