function [ z ] = rhsfunction( v )
%RHSFUNCTION Summary of this function goes here
%   Detailed explanation goes here
z=8*pi^2*sin(2*pi*v(:,1)).*cos(2*pi*v(:,2));

end

