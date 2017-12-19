syms a b s t Dst nu Dxr Dyr D1r Dxyr p;
NT=[-(s-1)*(t-1)*(s^2+s+t^2+t-2)/8;
    -b*(s-1)*(t-1)^2*(t+1);
    a*(s-1)^2*(s+1)*(t-1)/8;
    (s+1)*(t-1)*(s^2-s+t^2+t-2)/8;
    b*(s+1)*(t-1)^2*(t+1)/8;
    a*(s-1)*(s+1)^2*(t-1)/8;
    -(s+1)*(t+1)*(s^2-s+t^2-t-2)/8;
    b*(s+1)*(t-1)*(t+1)^2/8;
    -a*(s-1)*(s+1)^2*(t+1)/8;
    (s-1)*(t+1)*(s^2+s+t^2-t-2)/8;
    -b*(s-1)*(t-1)*(t+1)^2/8;
    -a*(s-1)^2*(s+1)*(t+1)/8];
Bnew=[-diff(NT,s,2)/a^2 -diff(NT,t,2)/b^2 -2*diff(diff(NT,s,1),t,1)/a/b];
%Cmat=Dst*[1 nu 0;nu 1 0; 0 0 (1-nu)/2];
Cmat=[Dxr D1r 0; D1r Dyr 0;0 0 Dxyr];
%Cmat=[1 nu 0;nu 1 0; 0 0 (1-nu)/2];
Koriginal=Bnew*Cmat*conj(Bnew');
%Kall=a*b*int(int(Bnew*Cmat*conj(Bnew'),t,-1,1),s,-1,1);
Kall=a*b*int(int(Koriginal,t,-1,1),s,-1,1);
check1=diff(Kall,Dxr)
check2=diff(Kall,Dyr)
check3=diff(Kall,D1r)
check4=diff(Kall,Dxyr)