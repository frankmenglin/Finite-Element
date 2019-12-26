function [ z ] = u_function( v )
%U_FUNCTION Summary of this function goes here
%   Detailed explanation goes here
z=sin(2*pi*v(:,1)).*cos(2*pi*v(:,2));

end

