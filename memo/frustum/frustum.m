function frustum
syms B BT eta C S l r z;

B(1,1:6) = [-1/l, -z/l/l*(12*eta - 6), -z/l*(6*eta - 4), 1/l, -z/l*(-12*eta + 6), -z/l*(6*eta - 2)];
B(2,1) = (1-eta)*S/r;
B(2,2) = 1/r*(C*(2*eta*eta*eta - 3*eta*eta + 1) - S/l*(6*eta*eta - 6*eta)); 
B(2,3) = l/r*(C*(eta*eta*eta - 2*eta*eta + eta) - S/l*(3*eta*eta - 4*eta + 1));
B(2,4) = eta/r*S;
B(2,5) = 1/r*(C*(- 2*eta*eta*eta + 3*eta*eta) - S/l*( - 6*eta*eta + 6*eta));
B(2,6) = l/r*(C*(eta*eta*eta - eta*eta) - S/l*(3*eta*eta - 2*eta));

for i = 1:2
    for j = 1:6
        BT(j,i) = B(i,j);
    end
end

syms k nu D;
D = k.*[1 nu
        nu 1];
syms Ke;
Ke = BT * D * B;

syms h T TT;
f2 = int(2*pi*l*int(r.*Ke,eta,0,1),z,-h/2,h/2);

T = [C S 0 0 0 0
     -S C 0 0 0 0 
     0 0 1 0 0 0 
     0 0 0 C S 0
     0 0 0 -S C 0 
     0 0 0 0 0 1];
for i = 1:6
    for j = 1:6
        TT(j,i) = T(i,j);
    end
end

ff = TT * f2 * T;

syms Matrix Matrixp0p;
Matrix = sym(zeros(1,21));
s = 0;
for j = 1:6
    for i = j:-1:1
        s = s + 1;
        Matrix(1,s) = ff(i,j);
    end
end

ccode(Matrix(1,:))
end