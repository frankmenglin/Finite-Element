function [ At,area ] = localstiffness( p )
%LOCALSTI Summary of this function goes here
%   Detailed explanation goes here
 At = zeros(3,3);
 B = [p(1,:)-p(3,:); p(2,:)-p(3,:)];
 G = [[1,0]',[0,1]',[-1,-1]'];
 area = 0.5*abs(det(B));
 for i = 1:3
 for j = 1:3
 At(i,j) = area*((B\G(:,i))'*(B\G(:,j)));
 end
 end

end

