
function[Y,x,C,D,BCD,B,Sh,Sv] = tireModel(alpha,gamma,Fz)

%Parameters
a0 = 1.30;
a1 = -7 * 10^-2;
a2 = 1.10;
a3 = 1.18; 
a4 = 7.80 ; 
a5 = 0.00; 
a6 = -0.20; 
a7 = 2.4 * 10^-2; 
a8 = 2.53 * 10^-2;
a9 = 0.00;
a10 = 0.00; 
a11 = 3.57 * 10^-3; 
a12 = 0.00;
a13 = 0.00; 

%Equations
Sh = (a8*gamma) + (a9*Fz) + a10; 
E = (a6*Fz) + a7; 
x = alpha + Sh; 
C = a0; 
D = (a1*Fz + a2)*Fz; 
BCD = a3*sin*(2*arctan(Fz/a4))*(1-a5*abs(gamma));
B = BCD/(C*D); 
Sv = (a11*Fz*gamma) + (a12*Fz) + a13; 
Y = D*sin(C*arctan*(B*x)-E*(B*x-arctan*(B*x)))+Sv; 

end 
