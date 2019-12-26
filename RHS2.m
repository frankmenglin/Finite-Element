function [ b ] = RHS2( node,elem )
%RIGHTHANDSIDE Summary of this function goes here
%   Input node,elem, function, compute the right hand side vector

 N = size(node,1);

 mid1 = (node(elem(:,2),:)+node(elem(:,3),:))/2;
 mid2 = (node(elem(:,3),:)+node(elem(:,1),:))/2;
 mid3 = (node(elem(:,1),:)+node(elem(:,2),:))/2; %Compute Midpoints
 
 %area =1.0/(size(elem,1)); %Compute Area of Triangles, only for unit square mesh
 
 bt1 = area.*(rhsfunction2(mid2)+rhsfunction2(mid3))/6;
 bt2 = area.*(rhsfunction2(mid3)+rhsfunction2(mid1))/6;
 bt3 = area.*(rhsfunction2(mid1)+rhsfunction2(mid2))/6;
 

 b = accumarray(elem(:),[bt1;bt2;bt3],[N 1]);

end


