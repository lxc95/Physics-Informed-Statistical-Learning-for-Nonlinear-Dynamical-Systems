function []=PiledData(r,c1_F_ROM_new,c1_StateData_displacement_ROM_new,c2_F_ROM_new,c2_StateData_displacement_ROM_new,c3_F_ROM_new,c3_StateData_displacement_ROM_new,c4_F_ROM_new,c4_StateData_displacement_ROM_new,c5_F_ROM_new,c5_StateData_displacement_ROM_new,c6_F_ROM_new,c6_StateData_displacement_ROM_new,c7_F_ROM_new,c7_StateData_displacement_ROM_new,c8_F_ROM_new,c8_StateData_displacement_ROM_new,c9_F_ROM_new,c9_StateData_displacement_ROM_new,c10_F_ROM_new,c10_StateData_displacement_ROM_new,c11_F_ROM_new,c11_StateData_displacement_ROM_new,c12_F_ROM_new,c12_StateData_displacement_ROM_new,c13_F_ROM_new,c13_StateData_displacement_ROM_new,c14_F_ROM_new,c14_StateData_displacement_ROM_new,c15_F_ROM_new,c15_StateData_displacement_ROM_new,c16_F_ROM_new,c16_StateData_displacement_ROM_new,c17_F_ROM_new,c17_StateData_displacement_ROM_new,c18_F_ROM_new,c18_StateData_displacement_ROM_new,c19_F_ROM_new,c19_StateData_displacement_ROM_new,c20_F_ROM_new,c20_StateData_displacement_ROM_new,c21_F_ROM_new,c21_StateData_displacement_ROM_new,c22_F_ROM_new,c22_StateData_displacement_ROM_new,c23_F_ROM_new,c23_StateData_displacement_ROM_new,c24_F_ROM_new,c24_StateData_displacement_ROM_new,c25_F_ROM_new,c25_StateData_displacement_ROM_new,c26_F_ROM_new,c26_StateData_displacement_ROM_new,c27_F_ROM_new,c27_StateData_displacement_ROM_new,c28_F_ROM_new,c28_StateData_displacement_ROM_new,c29_F_ROM_new,c29_StateData_displacement_ROM_new,c30_F_ROM_new,c30_StateData_displacement_ROM_new,c31_F_ROM_new,c31_StateData_displacement_ROM_new,c32_F_ROM_new,c32_StateData_displacement_ROM_new,c33_F_ROM_new,c33_StateData_displacement_ROM_new,c34_F_ROM_new,c34_StateData_displacement_ROM_new,c35_F_ROM_new,c35_StateData_displacement_ROM_new,c1_F_ROM_new_cum2,c2_F_ROM_new_cum2,c3_F_ROM_new_cum2,c4_F_ROM_new_cum2,c5_F_ROM_new_cum2,c6_F_ROM_new_cum2,c7_F_ROM_new_cum2,c8_F_ROM_new_cum2,c9_F_ROM_new_cum2,c10_F_ROM_new_cum2,c11_F_ROM_new_cum2,c12_F_ROM_new_cum2,c13_F_ROM_new_cum2,c14_F_ROM_new_cum2,c15_F_ROM_new_cum2,c16_F_ROM_new_cum2,c17_F_ROM_new_cum2,c18_F_ROM_new_cum2,c19_F_ROM_new_cum2,c20_F_ROM_new_cum2,c21_F_ROM_new_cum2,c22_F_ROM_new_cum2,c23_F_ROM_new_cum2,c24_F_ROM_new_cum2,c25_F_ROM_new_cum2,c26_F_ROM_new_cum2,c27_F_ROM_new_cum2,c28_F_ROM_new_cum2,c29_F_ROM_new_cum2,c30_F_ROM_new_cum2,c31_F_ROM_new_cum2,c32_F_ROM_new_cum2,c33_F_ROM_new_cum2,c34_F_ROM_new_cum2,c35_F_ROM_new_cum2)

num_POD = r;
for jk = 1:35
    dataset = [];
    dataset_name=strcat('c',num2str(jk),'_StateData_displacement_ROM_new');
    dataset = eval(dataset_name);
    temp = [];
    for mk = 1:num_POD
        temp = [temp dataset(mk,1:400)];
    end
    dataset_name2=strcat('c',num2str(jk),'_StateData_displacement_ROM_new_piled');
    assignin('base',dataset_name2,temp)
end

for jk = 1:35
    dataset = [];
    dataset_name=strcat('c',num2str(jk),'_F_ROM_new_cum2');
    dataset = eval(dataset_name);
    temp = [];
    for mk = 1:num_POD
        temp = [temp dataset(mk,1:400)];
    end
    dataset_name2=strcat('c',num2str(jk),'_F_ROM_new_cum2_piled');
    assignin('base',dataset_name2,temp)
end