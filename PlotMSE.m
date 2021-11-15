function []=PlotMSE(Uhat,full_predict_displacement_POD,c35_StateData_displacement)
full_predict_U = Uhat*full_predict_displacement_POD;
%MSE error all
figure
relative_MRE_error = zeros(1,200);
for i=1:200
    relative_MRE_error(1,i) = sum(abs(full_predict_U(:,i+200) - c35_StateData_displacement(:,i+200-30)))/sum(abs(c35_StateData_displacement(:,i+200-30)));
end
plot([4.02:0.02:8],relative_MRE_error*100,'-k','linewidth',2);
legend('MSE(%) of predicted snapshot-trajectory');
xlabel('Times (ms)','FontSize',18);
ylabel('MSE(%)','FontSize',18);
yLim = get(gca,'YLim');
set(gca,'YLim', [5,35]);
set(gca, 'FontName', 'Times New Roman')
set(gca,'FontSize',18)