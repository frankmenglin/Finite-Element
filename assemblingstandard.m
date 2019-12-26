function [ A ] = assemblingstandard(node,elem)
%ASSEMBLINGSTANDARD Summary of this function goes here
%   Detailed explanation goes here
N=size(node,1); NT=size(elem,1);
%A=zeros(N,N); 
A = sparse(N,N);
for t=1:NT
At=localstiffness(node(elem(t,:),:));
for i=1:3
for j=1:3
A(elem(t,i),elem(t,j))=A(elem(t,i),elem(t,j))+At(i,j);
end
end
end

end

