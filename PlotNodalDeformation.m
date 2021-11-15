function []=PlotNodalDeformation(gf,r,Uhat,full_predict_displacement_POD,c35_StateData_displacement,displacementFPCAbasis,covariance_FPCA_displacement)
full_predict_U = Uhat*full_predict_displacement_POD;
pre_full = full_predict_U;
nn=20;mm=12;% the truncation size in FPCA, POD 30
%plot the comparison of displacement
figure
%gf = 4999;
x_axis = [0.02:0.02:8];
plot(x_axis,pre_full(gf,:),'-r','linewidth',4);
hold on
plot(x_axis,[zeros(1,30) c35_StateData_displacement(gf,1:370)],'-.b','linewidth',4);
%calculate the variance
delta = zeros(1,400);
for t_step_predict = 1:400
    basis_special = displacementFPCAbasis(t_step_predict:400:(400*r),1:nn);
    covariance_POD_displacement = basis_special*covariance_FPCA_displacement*basis_special';
    variance_predict = zeros(40086,1);
    for i=1:40086
        variance_predict(i,1) = Uhat(i,:)*covariance_POD_displacement*Uhat(i,:)';
    end
    tem = sqrt(variance_predict);
    delta(t_step_predict) = tem(gf);
end
plot(x_axis,pre_full(gf,:)+1.96*delta,'--', 'Color', [.3 .3 .3],'linewidth',2);
plot(x_axis,pre_full(gf,:)-1.96*delta,'--', 'Color', [.3 .3 .3],'linewidth',2);
% xlim([-3 100])
% ylim([-20 6])
box off
h1 = legend('predicted','actual','95% confidence interval');
set(h1,'box','off');
xlabel('Time (ms)','FontSize',18);
ylabel('X-Displacement at node 1667','FontSize',18);
set(gca, 'FontName', 'Times New Roman')
set(gca,'FontSize',19)
