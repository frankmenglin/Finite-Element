function [ z ] = rhsfunction2sq( v )
%RHSFUNCTION2SQ Summary of this function goes here
%   Detailed explanation goes here
z=8*pi^2*cos(2*pi*v(:,1)).*cos(2*pi*v(:,2));

end

