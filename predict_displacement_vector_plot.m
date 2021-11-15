function [Mean_predict_displacement,Var_predict_displacement,temp_n] = predict_displacement_vector_plot(r,t_step_predict,nn,displacementFPCAbasis,covariance_FPCA_displacement,full_predict_displacement_POD,Uhat)
%t_step_predict = 400;
myfunc = @(x,y,z) sqrt(x.^2+y.^2+z.^2);
myfunc2 = @(x,y,z) x.^2+y.^2+z.^2;

basis_special = displacementFPCAbasis(t_step_predict:400:(400*r),1:nn);
covariance_POD_displacement = basis_special*covariance_FPCA_displacement*basis_special';
full_predict_U = Uhat*full_predict_displacement_POD(:,t_step_predict);
Mean_predict_displacement = zeros(13362,1);
Var_predict_displacement = zeros(13362,1);

temp_n = 0;

for j = 1:13362
    mu = [full_predict_U((j-1)*3+1) full_predict_U((j-1)*3+2) full_predict_U((j-1)*3+3)];
    
%     sigma = [ Uhat((j-1)*3+1,:)*covariance_POD_displacement*Uhat((j-1)*3+1,:)' Uhat((j-1)*3+1,:)*covariance_POD_displacement*Uhat((j-1)*3+2,:)' Uhat((j-1)*3+1,:)*covariance_POD_displacement*Uhat((j-1)*3+3,:)' ;
%         Uhat((j-1)*3+2,:)*covariance_POD_displacement*Uhat((j-1)*3+1,:)' Uhat((j-1)*3+2,:)*covariance_POD_displacement*Uhat((j-1)*3+2,:)' Uhat((j-1)*3+2,:)*covariance_POD_displacement*Uhat((j-1)*3+3,:)' ;
%         Uhat((j-1)*3+3,:)*covariance_POD_displacement*Uhat((j-1)*3+1,:)' Uhat((j-1)*3+3,:)*covariance_POD_displacement*Uhat((j-1)*3+2,:)' Uhat((j-1)*3+3,:)*covariance_POD_displacement*Uhat((j-1)*3+3,:)' ];
%     
    A_1 = Uhat((j-1)*3+1,:)*covariance_POD_displacement*Uhat((j-1)*3+1,:)';
    A_2 = Uhat((j-1)*3+1,:)*covariance_POD_displacement*Uhat((j-1)*3+2,:)';
    A_3 = Uhat((j-1)*3+1,:)*covariance_POD_displacement*Uhat((j-1)*3+3,:)';
    A_4 = Uhat((j-1)*3+2,:)*covariance_POD_displacement*Uhat((j-1)*3+2,:)';
    A_5 = Uhat((j-1)*3+2,:)*covariance_POD_displacement*Uhat((j-1)*3+3,:)';
    A_6 = Uhat((j-1)*3+3,:)*covariance_POD_displacement*Uhat((j-1)*3+3,:)';
    sigma = [A_1 A_2 A_3;
             A_2 A_4 A_5;
             A_3 A_5 A_6];
    tf = issymmetric(sigma);
     d = eig(sigma);
%d=0;
    issemidef = all(d > 0);
    if issemidef == 0 || tf == 0
        sigma = [ Uhat((j-1)*3+1,:)*covariance_POD_displacement*Uhat((j-1)*3+1,:)' 0 0 ;
            0 Uhat((j-1)*3+2,:)*covariance_POD_displacement*Uhat((j-1)*3+2,:)' 0 ;
            0 0 Uhat((j-1)*3+3,:)*covariance_POD_displacement*Uhat((j-1)*3+3,:)' ];
        % sigma = [ Uhat((j-1)*3+1,:)*covariance_POD_displacement*Uhat((j-1)*3+1,:)' Uhat((j-1)*3+2,:)*covariance_POD_displacement*Uhat((j-1)*3+1,:)' Uhat((j-1)*3+3,:)*covariance_POD_displacement*Uhat((j-1)*3+1,:)' ;
        %         Uhat((j-1)*3+2,:)*covariance_POD_displacement*Uhat((j-1)*3+1,:)' Uhat((j-1)*3+2,:)*covariance_POD_displacement*Uhat((j-1)*3+2,:)' Uhat((j-1)*3+3,:)*covariance_POD_displacement*Uhat((j-1)*3+2,:)' ;
        %         Uhat((j-1)*3+3,:)*covariance_POD_displacement*Uhat((j-1)*3+1,:)' Uhat((j-1)*3+3,:)*covariance_POD_displacement*Uhat((j-1)*3+2,:)' Uhat((j-1)*3+3,:)*covariance_POD_displacement*Uhat((j-1)*3+3,:)' ];
    temp_n = temp_n + 1;
    end
    random_x = mvnrnd(mu,sigma,20000)';
    temp_sum = myfunc(random_x(1,:),random_x(2,:),random_x(3,:));
    temp_sum_square = myfunc2(random_x(1,:),random_x(2,:),random_x(3,:));
    mean_temp_sum_square = mean(temp_sum_square);
    Mean_predict_displacement(j) = mean(temp_sum);
    Var_predict_displacement(j) = sqrt(mean_temp_sum_square-(mean(temp_sum))^2);
end