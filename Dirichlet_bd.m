function [ u ] = Dirichlet_bd(node,elem)
%DIRICHLET_BD Summary of this function goes here
%   Detailed explanation goes here
 A=assembling(node,elem);

 b=RightHandSide(node,elem);

 N = size(node,1);
 %isBdNode = false(N,1);
 Bdnode = findboundary(elem);
 FreeNode = setdiff([1:1:N]',Bdnode);
 %isBdNode(Dirichlet) = true;
 %bdNode = find(isBdNode);
 %freeNode = find(~isBdNode);
 u = zeros(N,1);
 u(Bdnode) = g_d_function(node(Bdnode,:));
 b = b - A*u;
 u(FreeNode)=A(FreeNode,FreeNode)\b(FreeNode);

end

