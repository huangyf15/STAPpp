function [Matrix0] = frustum
syms B BT eta;
syms C S l r0 k nu h;
syms xi1 xi2 weight1 weight2;
syms prmeter;

% prmeter = 2*pi*l
% k = E * h / (1 - nu*nu)

% C = sqrt(2)/2;
% S = sqrt(2)/2;
% r0 = 0.1;
% l = 0.2;
% xi1 = 0.8611363116;
% xi2 = 0.3399810436;
% weight1 = 0.3478548451;
% weight2 = 0.6521451549;
% k = 1000;
% nu = 0.2;
% h = 0.2;
% prmeter = 2*pi*l;

r = r0 + l * S * eta;
B(1,1) = -1/l;
B(1,4) = 1/l;
B(2,1) = (1 - eta) * S/r;
B(2,2) = (1 - 3*eta*eta + 2*eta*eta*eta) * C/r;
B(2,3) = l * (eta - 2*eta*eta + eta*eta*eta) * C/r;
B(2,4) = eta * S/r;
B(2,5) = (3*eta*eta - 2*eta*eta*eta) * C/r;
B(2,6) = l * (-eta*eta + eta*eta*eta) * C/r;
B(3,2) = -(-6 + 12*eta) / l / l;
B(3,3) = -(-4 + 6*eta) / l;
B(3,5) = -(6 - 12*eta) / l / l;
B(3,6) = -(-2 + 6*eta) / l;
B(4,2) = -(-6*eta + 6*eta*eta) * S/r/l;
B(4,3) = -(1 - 4*eta + 3*eta*eta) * S/r;
B(4,5) = -(6*eta - 6*eta*eta) * S/r/l;
B(4,6) = -(-2*eta + 3*eta*eta) * S/r;

for i = 1:4
    for j = 1:6
        BT(j,i) = B(i,j);
    end
end

% Calculate element stiffness matrix
D = sym(zeros(4));
D(1:2,1:2) = k.*[1 nu 
                 nu 1];
D(3:4,3:4) = h*h/12 * D(1:2,1:2);
Ke = BT * D * B;

% Numerical integration
f0 = sym(zeros(6));
for i = 1:6
    for j = 1:6
        f0(i,j) = weight1 * subs(r*Ke(i,j),eta,(1+xi1)*0.5) ...
                  + weight1 * subs(r*Ke(i,j),eta,(1-xi1)*0.5) ...
                  + weight2 * subs(r*Ke(i,j),eta,(1+xi2)*0.5) ...
                  + weight2 * subs(r*Ke(i,j),eta,(1-xi2)*0.5);
    end
end

% Transform to the background coordinate
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
ff = prmeter*TT*f0*T;

% Transform to the symbolic form
syms Matrix Matrixp0p;
Matrix = sym(zeros(1,21));
s = 0;
for j = 1:6
    for i = j:-1:1
        s = s + 1;
        Matrix(1,s) = ff(i,j);
    end
end

Matrix0 = ccode(Matrix(1,:));
end