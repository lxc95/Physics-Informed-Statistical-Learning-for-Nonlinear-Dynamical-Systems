function []=PlotPODbasis(temp)
load('fem_node_element.mat');load('Coordinate_XYZ.mat');
%temp = Uhat(:,10)';
Displacement_X_mm=[temp(:,1:3:40086)];
Displacement_Y_mm=[temp(:,2:3:40086)];
Displacement_Z_mm=[temp(:,3:3:40086)];

Real_coordinate_x=zeros(1,13362);
Real_coordinate_y=zeros(1,13362);
Real_coordinate_z=zeros(1,13362);
Real_coordinate_x=Displacement_X_mm+Coordinate_X_mm(1,:);
Real_coordinate_y=Displacement_Y_mm+Coordinate_Y_mm(1,:);
Real_coordinate_z=Displacement_Z_mm+Coordinate_Z_mm(1,:);
Displacement_all=zeros(1,13362);
for i=1:13362
    Displacement_all(1,i)=sqrt((Displacement_X_mm(i))^2+(Displacement_Y_mm(i))^2+(Displacement_Z_mm(i))^2);
end
%% plot of Displacement
%contour of Displacement
%1,component
component=[];
%figure;
%extended Real_coordinate_all
extended_Real_Displacement_all=[];
extended_Real_Displacement_all=zeros(2545746,1);
extended_Real_Displacement_all=[extended_Real_Displacement_all;Displacement_all'];
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
axis off ;
colormap(jet(15));
view(40,50);
%set(gcf,'Units','centimeters','position',[5 5 100 100]);
set(gcf,'Units','centimeters');

h = colorbar;
c = jet(15);
colormap(c);
set(get(h,'label'),'string','Displacement(mm)');
set(gca, 'FontName', 'Times New Roman')
set(gca,'FontSize',48)
