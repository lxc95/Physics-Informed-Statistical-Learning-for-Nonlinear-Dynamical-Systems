function [Uhat,Sig,c1_F_ROM,c1_StateData_displacement_ROM,c2_F_ROM,c2_StateData_displacement_ROM,c3_F_ROM,c3_StateData_displacement_ROM,c4_F_ROM,c4_StateData_displacement_ROM,c5_F_ROM,c5_StateData_displacement_ROM,c6_F_ROM,c6_StateData_displacement_ROM,c7_F_ROM,c7_StateData_displacement_ROM,c8_F_ROM,c8_StateData_displacement_ROM,c9_F_ROM,c9_StateData_displacement_ROM,c10_F_ROM,c10_StateData_displacement_ROM,c11_F_ROM,c11_StateData_displacement_ROM,c12_F_ROM,c12_StateData_displacement_ROM,c13_F_ROM,c13_StateData_displacement_ROM,c14_F_ROM,c14_StateData_displacement_ROM,c15_F_ROM,c15_StateData_displacement_ROM,c16_F_ROM,c16_StateData_displacement_ROM,c17_F_ROM,c17_StateData_displacement_ROM,c18_F_ROM,c18_StateData_displacement_ROM,c19_F_ROM,c19_StateData_displacement_ROM,c20_F_ROM,c20_StateData_displacement_ROM,c21_F_ROM,c21_StateData_displacement_ROM,c22_F_ROM,c22_StateData_displacement_ROM,c23_F_ROM,c23_StateData_displacement_ROM,c24_F_ROM,c24_StateData_displacement_ROM,c25_F_ROM,c25_StateData_displacement_ROM,c26_F_ROM,c26_StateData_displacement_ROM,c27_F_ROM,c27_StateData_displacement_ROM,c28_F_ROM,c28_StateData_displacement_ROM,c29_F_ROM,c29_StateData_displacement_ROM,c30_F_ROM,c30_StateData_displacement_ROM,c31_F_ROM,c31_StateData_displacement_ROM,c32_F_ROM,c32_StateData_displacement_ROM,c33_F_ROM,c33_StateData_displacement_ROM,c34_F_ROM,c34_StateData_displacement_ROM,c35_F_ROM,c35_StateData_displacement_ROM]=POD(r,full_M_all_lumped)
load('fem_node_element.mat');load('alldata_Force_Displacement.mat');
%reconstruct element matrix
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
% SVD
allStateData_dis = [c1_StateData_displacement,c2_StateData_displacement,c3_StateData_displacement,c4_StateData_displacement,c5_StateData_displacement,c6_StateData_displacement,c7_StateData_displacement,c8_StateData_displacement,c9_StateData_displacement,c10_StateData_displacement,c11_StateData_displacement,c12_StateData_displacement,c13_StateData_displacement,c14_StateData_displacement,c15_StateData_displacement,c16_StateData_displacement,c17_StateData_displacement,c18_StateData_displacement,c19_StateData_displacement,c20_StateData_displacement,c21_StateData_displacement,c22_StateData_displacement,c23_StateData_displacement,c24_StateData_displacement,c25_StateData_displacement,c26_StateData_displacement,c27_StateData_displacement,c28_StateData_displacement,c29_StateData_displacement,c30_StateData_displacement,c31_StateData_displacement,c32_StateData_displacement,c33_StateData_displacement,c34_StateData_displacement,c35_StateData_displacement];
Xp = allStateData_dis(:,1:3:end);%because the size of column is too large
Xp_new = sqrt(full_M_all_lumped)*Xp;
[U,Sig,V] = svd(Xp_new,'econ');
Uhat = inv(sqrt(full_M_all_lumped))*U(:,1:r);
%.................................................................................................................................
%projection with POD basis
InputData = c1_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c1_F = F_1 + F_2;
c1_F_ROM = Uhat'*c1_F;
c1_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c1_StateData_displacement;

InputData = c2_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c2_F = F_1 + F_2;
c2_F_ROM = Uhat'*c2_F;
c2_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c2_StateData_displacement;

InputData = c3_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c3_F = F_1 + F_2;
c3_F_ROM = Uhat'*c3_F;
c3_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c3_StateData_displacement;

InputData = c4_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c4_F = F_1 + F_2;
c4_F_ROM = Uhat'*c4_F;
c4_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c4_StateData_displacement;

InputData = c5_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c5_F = F_1 + F_2;
c5_F_ROM = Uhat'*c5_F;
c5_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c5_StateData_displacement;

InputData = c6_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c6_F = F_1 + F_2;
c6_F_ROM = Uhat'*c6_F;
c6_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c6_StateData_displacement;

InputData = c7_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c7_F = F_1 + F_2;
c7_F_ROM = Uhat'*c7_F;
c7_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c7_StateData_displacement;

InputData = c8_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c8_F = F_1 + F_2;
c8_F_ROM = Uhat'*c8_F;
c8_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c8_StateData_displacement;

InputData = c9_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c9_F = F_1 + F_2;
c9_F_ROM = Uhat'*c9_F;
c9_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c9_StateData_displacement;

InputData = c10_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c10_F = F_1 + F_2;
c10_F_ROM = Uhat'*c10_F;
c10_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c10_StateData_displacement;

InputData = c11_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c11_F = F_1 + F_2;
c11_F_ROM = Uhat'*c11_F;
c11_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c11_StateData_displacement;

InputData = c12_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c12_F = F_1 + F_2;
c12_F_ROM = Uhat'*c12_F;
c12_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c12_StateData_displacement;

InputData = c13_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c13_F = F_1 + F_2;
c13_F_ROM = Uhat'*c13_F;
c13_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c13_StateData_displacement;

InputData = c14_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c14_F = F_1 + F_2;
c14_F_ROM = Uhat'*c14_F;
c14_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c14_StateData_displacement;

InputData = c15_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c15_F = F_1 + F_2;
c15_F_ROM = Uhat'*c15_F;
c15_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c15_StateData_displacement;

InputData = c16_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c16_F = F_1 + F_2;
c16_F_ROM = Uhat'*c16_F;
c16_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c16_StateData_displacement;

InputData = c17_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c17_F = F_1 + F_2;
c17_F_ROM = Uhat'*c17_F;
c17_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c17_StateData_displacement;

InputData = c18_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c18_F = F_1 + F_2;
c18_F_ROM = Uhat'*c18_F;
c18_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c18_StateData_displacement;

InputData = c19_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c19_F = F_1 + F_2;
c19_F_ROM = Uhat'*c19_F;
c19_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c19_StateData_displacement;

InputData = c20_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c20_F = F_1 + F_2;
c20_F_ROM = Uhat'*c20_F;
c20_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c20_StateData_displacement;

InputData = c21_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c21_F = F_1 + F_2;
c21_F_ROM = Uhat'*c21_F;
c21_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c21_StateData_displacement;

InputData = c22_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c22_F = F_1 + F_2;
c22_F_ROM = Uhat'*c22_F;
c22_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c22_StateData_displacement;

InputData = c23_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c23_F = F_1 + F_2;
c23_F_ROM = Uhat'*c23_F;
c23_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c23_StateData_displacement;

InputData = c24_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c24_F = F_1 + F_2;
c24_F_ROM = Uhat'*c24_F;
c24_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c24_StateData_displacement;

InputData = c25_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c25_F = F_1 + F_2;
c25_F_ROM = Uhat'*c25_F;
c25_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c25_StateData_displacement;

InputData = c26_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c26_F = F_1 + F_2;
c26_F_ROM = Uhat'*c26_F;
c26_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c26_StateData_displacement;

InputData = c27_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c27_F = F_1 + F_2;
c27_F_ROM = Uhat'*c27_F;
c27_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c27_StateData_displacement;

InputData = c28_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c28_F = F_1 + F_2;
c28_F_ROM = Uhat'*c28_F;
c28_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c28_StateData_displacement;

InputData = c29_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c29_F = F_1 + F_2;
c29_F_ROM = Uhat'*c29_F;
c29_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c29_StateData_displacement;

InputData = c30_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c30_F = F_1 + F_2;
c30_F_ROM = Uhat'*c30_F;
c30_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c30_StateData_displacement;

InputData = c31_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c31_F = F_1 + F_2;
c31_F_ROM = Uhat'*c31_F;
c31_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c31_StateData_displacement;

InputData = c32_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c32_F = F_1 + F_2;
c32_F_ROM = Uhat'*c32_F;
c32_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c32_StateData_displacement;

InputData = c33_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c33_F = F_1 + F_2;
c33_F_ROM = Uhat'*c33_F;
c33_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c33_StateData_displacement;

InputData = c34_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c34_F = F_1 + F_2;
c34_F_ROM = Uhat'*c34_F;
c34_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c34_StateData_displacement;

InputData = c35_InputData;
F_1 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell3P',13362,InputData);
F_2 = ForceAssembler3Dshell(num_and_coordinate(:,2:4)',newshell4P',13362,InputData);
c35_F = F_1 + F_2;
c35_F_ROM = Uhat'*c35_F;
c35_StateData_displacement_ROM = Uhat'*full_M_all_lumped*c35_StateData_displacement;


