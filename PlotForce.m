function []=PlotForce(c1_F_ROM_new_cum2,c2_F_ROM_new_cum2,c3_F_ROM_new_cum2,c4_F_ROM_new_cum2,c5_F_ROM_new_cum2,c6_F_ROM_new_cum2,c7_F_ROM_new_cum2,c8_F_ROM_new_cum2,c9_F_ROM_new_cum2,c10_F_ROM_new_cum2,c11_F_ROM_new_cum2,c12_F_ROM_new_cum2,c13_F_ROM_new_cum2,c14_F_ROM_new_cum2,c15_F_ROM_new_cum2,c16_F_ROM_new_cum2,c17_F_ROM_new_cum2,c18_F_ROM_new_cum2,c19_F_ROM_new_cum2,c20_F_ROM_new_cum2,c21_F_ROM_new_cum2,c22_F_ROM_new_cum2,c23_F_ROM_new_cum2,c24_F_ROM_new_cum2,c25_F_ROM_new_cum2,c26_F_ROM_new_cum2,c27_F_ROM_new_cum2,c28_F_ROM_new_cum2,c29_F_ROM_new_cum2,c30_F_ROM_new_cum2,c31_F_ROM_new_cum2,c32_F_ROM_new_cum2,c33_F_ROM_new_cum2,c34_F_ROM_new_cum2,c35_F_ROM_new_cum2)
figure
for jk = 1:35
    dataset = [];
    dataset_name=strcat('c',num2str(jk),'_F_ROM_new_cum2');
    dataset = eval(dataset_name);
    tempset = dataset(:,1:400);
    
    
    t = (0:0.02:7.98)';
    tMat = repmat(t, 1, 10); %// For plot3
    tempset = tempset';
    
    index = 1:1:10;
    indexMat = repmat(index, numel(t), 1); %//For plot3
    plot3(tMat, indexMat, tempset(:,1:10)); %// Make all traces blue
    grid;
    xlabel('Time (ms)'); ylabel('POD index'); zlabel('POD scores of cumulative external force');
    view(40,40); %// Adjust viewing angle so you can clearly see data
    hold on
end
set(gca, 'FontName', 'Times New Roman')
set(gca,'FontSize',16)