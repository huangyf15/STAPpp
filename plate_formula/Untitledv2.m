syms a b s t Dst nu Dxr Dyr D1r Dxyr p;
NT=[-(s-1)*(t-1)*(s^2+s+t^2+t-2)/8;
    -(s-1)*(t-1)^2*(t+1);
    (s-1)^2*(s+1)*(t-1)/8;
    (s+1)*(t-1)*(s^2-s+t^2+t-2)/8;
    (s+1)*(t-1)^2*(t+1)/8;
    (s-1)*(s+1)^2*(t-1)/8;
    -(s+1)*(t+1)*(s^2-s+t^2-t-2)/8;
    (s+1)*(t-1)*(t+1)^2/8;
    -(s-1)*(s+1)^2*(t+1)/8;
    (s-1)*(t+1)*(s^2+s+t^2-t-2)/8;
    -(s-1)*(t-1)*(t+1)^2/8;
    -(s-1)^2*(s+1)*(t+1)/8];
Bnew=[-diff(NT,s,2) -diff(NT,t,2) -diff(diff(NT,s,1),t,1)];

syms xpsi xeta ypsi yeta psix psiy etax etay
J = [xpsi ypsi;xeta yeta];
invJ = inv(J);
%chg_mat = [invJ(1,1)^2 invJ(1,2)^2 2*invJ(1,1)*invJ(1,2);
%           invJ(2,1)^2 invJ(2,2)^2 2*invJ(2,1)*invJ(2,2);
%           2*invJ(1,1)*invJ(2,1) 2*invJ(1,2)*invJ(2,2) 2*(invJ(1,1)*invJ(2,2)+invJ(1,2)+invJ(2,1))];
chg_mat = [psix*psix etax*etax 2*psix*etax;
           psiy*psiy etay*etay 2*psiy*etay;
           2*psix*psiy 2*etax*etay 2*(psix*etay+psiy*etax)];


%Cmat=Dst*[1 nu 0;nu 1 0; 0 0 (1-nu)/2];
%Cmat=[Dxr D1r 0; D1r Dyr 0;0 0 Dxyr];
Cmat=[1 nu 0;nu 1 0; 0 0 (1-nu)/2];
Koriginal=Bnew*conj(chg_mat')*Cmat*chg_mat*conj(Bnew');
%Kall=a*b*int(int(Bnew*Cmat*conj(Bnew'),t,-1,1),s,-1,1);
Kall=a*b*int(int(Koriginal,t,-1,1),s,-1,1);

fileused=fopen('newint.txt','wt');
for j=1:12
    for i=j:-1:1
        fprintf(fileused,'%s%i%s%s%s\n','Matrix[',j*(j+1)/2-i,'] = ',char(Kall(i,j)),';');
    end
end
fclose(fileused);
        