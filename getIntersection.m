function xy = getIntersection(L1,L2)


x = (L1(2)*L2(3)-L2(2)*L1(3))/(L1(1)*L2(2)-L2(1)*L1(2)); 
y = (L1(3)*L2(1)-L2(3)*L1(1))/(L1(1)*L2(2)-L2(1)*L1(2)); 
xy = [-x -y];

%x = (b1*c2-b2*c1)/(a1*b2-a2*b1); 
%y = (c1*a2-c2*a1)/(a1*b2-a2*b1); 