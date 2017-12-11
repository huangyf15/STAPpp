syms xpsi xeta ypsi yeta
J = [xpsi ypsi;xeta yeta];
invJ = inv(J);
chg_mat = [invJ(1,1)^2 invJ(1,2)^2 2*invJ(1,1)*invJ(1,2);
           invJ(2,1)^2 invJ(2,2)^2 2*invJ(2,1)*invJ(2,2);
           invJ(1,1)*invJ(2,1) invJ(1,2)*invJ(2,2) invJ(1,1)*invJ(2,2)+invJ(1,2)+invJ(2,1)];
ck = subs(subs(subs(subs(chg_mat,xpsi,1),ypsi,0),xeta,0),yeta,1)