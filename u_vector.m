function [ u ] = u_vector( node,elem )
%U_VECTOR Summary of this function goes here
%   Detailed explanation goes here
N = size(node,1);

v1 = node(elem(:,1),:);
v2 = node(elem(:,2),:);
v3 = node(elem(:,3),:);

ut1 = u_function(v1)/3;
ut2 = u_function(v2)/3;
ut3 = u_function(v3)/3;

u = accumarray(elem(:),[ut1;ut2;ut3],[N 1]);

end

