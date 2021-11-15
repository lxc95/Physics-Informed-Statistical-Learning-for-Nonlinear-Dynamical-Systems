function []=PlotXstandarddeviation(Uhat,covariance_POD_displacement)
load('fem_node_element.mat');load('Coordinate_XYZ.mat');
variance_predict = zeros(40086,1);
for i=1:40086
    variance_predict(i,1) = Uhat(i,:)*covariance_POD_displacement*Uhat(i,:)';
end
% figure
% plot(sqrt(variance_predict));
%pre_full = sqrt(variance_predict)*1.96;
pre_full = sqrt(variance_predict);
temp = pre_full';
Displacement_X_mm=[temp(:,1:3:40086)];
Displacement_Y_mm=[temp(:,2:3:40086)];
Displacement_Z_mm=[temp(:,3:3:40086)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
for bn = 1
    num_step = bn;
    Real_coordinate_x=zeros(1,13362);
    Real_coordinate_y=zeros(1,13362);
    Real_coordinate_z=zeros(1,13362);
    Real_coordinate_x=Displacement_X_mm(num_step,:)+Coordinate_X_mm(1,:);
    Real_coordinate_y=Displacement_Y_mm(num_step,:)+Coordinate_Y_mm(1,:);
    Real_coordinate_z=Displacement_Z_mm(num_step,:)+Coordinate_Z_mm(1,:);
    Displacement_all=zeros(1,13362);
    for i=1:13362
        Displacement_all(1,i)=sqrt((Displacement_X_mm(num_step,i))^2+(Displacement_Y_mm(num_step,i))^2+(Displacement_Z_mm(num_step,i))^2);
    end
    %% plot of Displacement
    %contour of Displacement
    %1,component
    component=[];
    %figure;
    %extended Real_coordinate_all
    extended_Real_Displacement_all=[];
    extended_Real_Displacement_all=zeros(2545746,1);
    extended_Real_Displacement_all=[extended_Real_Displacement_all;abs(Displacement_X_mm(num_step,:)')];%the contour
    component = extended_Real_Displacement_all;

    %2,coordinates
    extendedShell=[];
    extendedShell=zeros(2545746,3);
    extendedShell=[extendedShell;[Real_coordinate_x' Real_coordinate_y' Real_coordinate_z']];
    coordinates=extendedShell;
    %3,nodes
    nodes=allshell(67:13033,3:6);
    %use plot function
    PlotFieldonMesh(coordinates,nodes,component) ; % Plot the component profile on mesh
    hold on
    f=allshell(1:66,3:5);
    coordinates=extendedShell;
    nodes=f;
    nnel=3;
    nel=66;
    X = zeros(nnel,nel) ;
    Y = zeros(nnel,nel) ;
    Z = zeros(nnel,nel) ;
    profile = zeros(nnel,nel) ;
    for iel=1:nel
        nd=nodes(iel,:);         % extract connected node for (iel)-th element
        X(:,iel)=coordinates(nd,1);    % extract x value of the node
        Y(:,iel)=coordinates(nd,2);    % extract y value of the node
        Z(:,iel)=coordinates(nd,3) ;   % extract z value of the node
        profile(:,iel) = component(nd') ; % extract component value of the node
    end
    % Plotting the FEM mesh and profile of the given component
    %figure
    fill3(X,Y,Z,profile)
    rotate3d on ;
    %title('Profile of component on Mesh') ;
    %axes('position', [0, 0, 0, 0.5]);
    axis off ;
    colormap(jet(15));
    view(40,50);
    scrsz = get(0,'ScreenSize');
    %set(gcf,'Units','centimeters','position',scrsz);
    set(gcf,'Units','centimeters');
%     p=1;
%     fig(p)=figure(p);
%     picturename=strcat('real_FOM_step_',num2str(num_step),'.jpg');
%     saveas(fig(p),picturename,'jpg');
    h = colorbar;
    c = jet(15);
    %colormap(c),caxis([0 28]);
    set(get(h,'label'),'string','X-Displacement(mm)');
    set(gca, 'FontName', 'Times New Roman')
    set(gca,'FontSize',36)
    %close all;
end

