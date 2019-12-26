%  Example: compute L2 error of piecewise linear interpolation

     [node,elem] = squaremesh([0,1,0,1],0.25);
     for k = 1:4
         exactu = inline('sin(pi*pxy(:,1)).*sin(pi*pxy(:,2))','pxy');
         uI = exactu(node);
         N(k) = size(node,1);
         err(k) = getL2error(node,elem,exactu,uI);
         [node,elem] = uniformrefine(node,elem);
     end
     showrate(N,err);