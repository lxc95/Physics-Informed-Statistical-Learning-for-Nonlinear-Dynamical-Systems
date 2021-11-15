function M = MassAssembler3Dshell(p,t,n)
%p is the node matrix; t is the element matrix; n is the dimension of M
if size(t,1)==3
    %assemble the triangle elements
    %np = size(p,2);
    nt = size(t,2);
    M = sparse(n,n);
    for K = 1:nt
        loc2glb = t(1:3,K);
        x = p(1,loc2glb);
        y = p(2,loc2glb);
        z = p(3,loc2glb);
        A = [x(1) y(1) z(1)];
        B = [x(2) y(2) z(2)];
        C = [x(3) y(3) z(3)];
        area = Area3D3P(A,B,C);
        MK1 = [2 1 1;
            1 2 1;
            1 1 2]/12*area;
        M(loc2glb,loc2glb) = M(loc2glb,loc2glb) + MK1;
    end
else
    %assemble the quadrilateral
    %np = size(p,2);
    nt = size(t,2);
    M = sparse(n,n);
    for K = 1:nt
        loc2glb = t(1:4,K);
        x = p(1,loc2glb);
        y = p(2,loc2glb);
        z = p(3,loc2glb);
        A = [x(1) y(1) z(1)];
        B = [x(2) y(2) z(2)];
        C = [x(3) y(3) z(3)];
        D = [x(4) y(4) z(4)];
        area = Area3D4P(A,B,C,D);
%         MK2 = [2 1 1 1;
%             1 2 1 1;
%             1 1 2 1;
%             1 1 1 2]/20*area;
        MK2 = [4 2 1 2;
               2 4 2 1;
               1 2 4 2;
               2 1 2 4]/36*area;
        M(loc2glb,loc2glb) = M(loc2glb,loc2glb) + MK2;
        %area
    end
    
end