 clear;
 exactu = inline('pxy(:,1).^2+pxy(:,2).^2-0.5','pxy');
 Du = inline('[2*pxy(:,1) 2*pxy(:,2)]','pxy');
 
 for k=1:4
 [node,elem] = circlemesh(0,0,1,0.4/2^k);
 A= assembling(node,elem);
 N(k) = size(node,1);
 mid1 = (node(elem(:,2),:)+node(elem(:,3),:))/2;
 mid2 = (node(elem(:,3),:)+node(elem(:,1),:))/2;
 mid3 = (node(elem(:,1),:)+node(elem(:,2),:))/2;
 
 e1 = (node(elem(:,2),:)-node(elem(:,1),:));
 e1 = [e1 zeros(size(e1,1),1)];
 e2 = (node(elem(:,3),:)-node(elem(:,1),:));
 e2 = [e2 zeros(size(e2,1),1)];
 area = abs((cross(e1,e2)))*1/2;
 area=area(:,3);

 
 
 bt1 = area.*(rhsfunction2(mid2)+rhsfunction2(mid3))/6;
 bt2 = area.*(rhsfunction2(mid3)+rhsfunction2(mid1))/6;
 bt3 = area.*(rhsfunction2(mid1)+rhsfunction2(mid2))/6;
 
 b = accumarray(elem(:),[bt1;bt2;bt3],[N(k) 1]);
 
 [bdNode,bdEdge,isBdNode,isBdElem] = findboundary(elem);
 [bdNode,bdEdge,isBdNode,isBdElem] = findboundary(elem);
 Nve = node(bdEdge(:,1),:) - node(bdEdge(:,2),:);
 edgeLength = sqrt(sum(Nve.^2,2));
 mid = (node(bdEdge(:,1),:) + node(bdEdge(:,2),:))/2;
 b = b + accumarray([bdEdge(:),ones(2*size(bdEdge,1),1)],repmat(edgeLength.*g_N(mid)/2,2,1),[N(k),1]);
 b = b - mean(b);
 u_h=A\b;
 u_h=u_h-mean(u_h);
 error_L2(k) = getL2error(node,elem,exactu,u_h);
 error_H1(k) = getH1error(node,elem,Du,u_h);

 end
 %%
 showrate(N,error_L2);
 %%
 showrate(N,error_H1);
 %%
  showresult(node,elem,u_h-mean(u_h));
