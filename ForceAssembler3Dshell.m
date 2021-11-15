function F = ForceAssembler3Dshell(p,t,n,input)
%p is the node matrix; t is the element matrix; n is the dimension of M;
%input is the force term
ndof = 3*size(p,2);% total number of degrees of freedom
F = zeros(ndof,size(input,2));
if size(t,1)==3
    dofs = zeros(9,1);%allocate element degrees of freedom
    %assemble the triangle elements
    %np = size(p,2);
    nt = size(t,2);
    %M = sparse(n,n);
    for time = 1:size(input,2)
        for K = 1:nt
            loc2glb = t(1:3,K);
            x = p(1,loc2glb);
            y = p(2,loc2glb);
            z = p(3,loc2glb);
            dofs(3:3:end) = 3*loc2glb;%element degrees of freedom
            dofs(2:3:end) = 3*loc2glb-1;
            dofs(1:3:end) = 3*loc2glb-2;
            A = [x(1) y(1) z(1)];
            B = [x(2) y(2) z(2)];
            C = [x(3) y(3) z(3)];
            area = Area3D3P(A,B,C);
            bk = input(dofs,time)/3*area;
            F(dofs,time) = F(dofs,time) + bk;
        end
    end
else
    dofs = zeros(12,1);%allocate element degrees of freedom
    %assemble the quadrilateral
    %np = size(p,2);
    nt = size(t,2);
    %M = sparse(n,n);
    for time = 1:size(input,2)
        for K = 1:nt
            loc2glb = t(1:4,K);
            x = p(1,loc2glb);
            y = p(2,loc2glb);
            z = p(3,loc2glb);
            dofs(3:3:end) = 3*loc2glb;%element degrees of freedom
            dofs(2:3:end) = 3*loc2glb-1;
            dofs(1:3:end) = 3*loc2glb-2;
            A = [x(1) y(1) z(1)];
            B = [x(2) y(2) z(2)];
            C = [x(3) y(3) z(3)];
            D = [x(4) y(4) z(4)];
            area = Area3D4P(A,B,C,D);
            bk = input(dofs,time)/4*area;
            F(dofs,time) = F(dofs,time) + bk;
        end
    end
end
% function F = ForceAssembler3Dshell(p,t,n,input)
% %p is the node matrix; t is the element matrix; n is the dimension of M;
% %input is the force term
% ndof = 3*size(p,2);% total number of degrees of freedom
% M = sparse(ndof,ndof);
% F = zeros(ndof,size(input,2));
% if size(t,1)==3
%     dofs = zeros(9,1);%allocate element degrees of freedom
%     %assemble the triangle elements
%     %np = size(p,2);
%     nt = size(t,2);
%     %M = sparse(n,n);
%     for K = 1:nt
%         loc2glb = t(1:3,K);
%         x = p(1,loc2glb);
%         y = p(2,loc2glb);
%         z = p(3,loc2glb);
%         dofs(3:3:end) = 3*loc2glb;%element degrees of freedom
%         dofs(2:3:end) = 3*loc2glb-1;
%         dofs(1:3:end) = 3*loc2glb-2;
%         A = [x(1) y(1) z(1)];
%         B = [x(2) y(2) z(2)];
%         C = [x(3) y(3) z(3)];
%         area = Area3D3P(A,B,C);
%         MK1 = [2 0 0 1 0 0 1 0 0;
%                0 2 0 0 1 0 0 1 0;
%                0 0 2 0 0 1 0 0 1;
%                1 0 0 2 0 0 1 0 0;
%                0 1 0 0 2 0 0 1 0;
%                0 0 1 0 0 2 0 0 1;
%                1 0 0 1 0 0 2 0 0;
%                0 1 0 0 1 0 0 2 0;
%                0 0 1 0 0 1 0 0 2]/12*area;
%         fK = input(dofs,:);
%         FK = MK1*fK;
%         M(dofs,dofs) = M(dofs,dofs) + MK1;
%         F(dofs,:) = F(dofs,:) + FK;
%     end
% else
%     dofs = zeros(12,1);%allocate element degrees of freedom
%     %assemble the quadrilateral
%     %np = size(p,2);
%     nt = size(t,2);
%     %M = sparse(n,n);
%     for K = 1:nt
%         loc2glb = t(1:4,K);
%         x = p(1,loc2glb);
%         y = p(2,loc2glb);
%         z = p(3,loc2glb);
%         dofs(3:3:end) = 3*loc2glb;%element degrees of freedom
%         dofs(2:3:end) = 3*loc2glb-1;
%         dofs(1:3:end) = 3*loc2glb-2;
%         A = [x(1) y(1) z(1)];
%         B = [x(2) y(2) z(2)];
%         C = [x(3) y(3) z(3)];
%         D = [x(4) y(4) z(4)];
%         area = Area3D4P(A,B,C,D);
% %         MK2 = [2 0 0 1 0 0 1 0 0 1 0 0;
% %                0 2 0 0 1 0 0 1 0 0 1 0;
% %                0 0 2 0 0 1 0 0 1 0 0 1;
% %                1 0 0 2 0 0 1 0 0 1 0 0;
% %                0 1 0 0 2 0 0 1 0 0 1 0;
% %                0 0 1 0 0 2 0 0 1 0 0 1;
% %                1 0 0 1 0 0 2 0 0 1 0 0;
% %                0 1 0 0 1 0 0 2 0 0 1 0;
% %                0 0 1 0 0 1 0 0 2 0 0 1;
% %                1 0 0 1 0 0 1 0 0 2 0 0;
% %                0 1 0 0 1 0 0 1 0 0 2 0;
% %                0 0 1 0 0 1 0 0 1 0 0 2]/20*area;
%         MK2 = [4 0 0 2 0 0 1 0 0 2 0 0;
%                0 4 0 0 2 0 0 1 0 0 2 0;
%                0 0 4 0 0 2 0 0 1 0 0 2;
%                2 0 0 4 0 0 2 0 0 1 0 0;
%                0 2 0 0 4 0 0 2 0 0 1 0;
%                0 0 2 0 0 4 0 0 2 0 0 1;
%                1 0 0 2 0 0 4 0 0 2 0 0;
%                0 1 0 0 2 0 0 4 0 0 2 0;
%                0 0 1 0 0 2 0 0 4 0 0 2;
%                2 0 0 1 0 0 2 0 0 4 0 0;
%                0 2 0 0 1 0 0 2 0 0 4 0;
%                0 0 2 0 0 1 0 0 2 0 0 4]/36*area;
%         fK = input(dofs,:);
%         FK = MK2*fK;
%         M(dofs,dofs) = M(dofs,dofs) + MK2;
%         F(dofs,:) = F(dofs,:) + FK;
%     end
%     
% end
