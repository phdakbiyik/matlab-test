clc;
clear all;
syms s
r = 1/sqrt(3)*[2*sin(s),(1/2)*cos(2*s),(1/2)*sin(2*s),(1/3)*cos(3*s), ...
                2*cos(s),(1/4)*sin(4*s),(1/4)*cos(4*s),(1/3)*sin(3*s)];
[G, dGds, A, C] = G2associatedframe(r, pi/4);

function [G, dGds, A, C] = G2associatedframe(r, s_val)
% Computes the G2 associated frame and its derivative 
% for a given unit speed split-octonion curve
%
% INPUTS:
%   r : a unit speed split octonionic curve r(s) 
%   Example r = [r1(s), r2(s), ..., r8(s)];
%   s_val  : evaluation point for numerical results (e.g., pi/4)
%
% OUTPUTS:
%   G      : 8x8 frame matrix at s_val
%   dGds   : derivative matrix B*G at s_val
%   A      :  A matrix at s_val
%   C      :  C matrix at s_val

% abbreviations: ei=epsilon_i, mi=mu_i, kai=kappa_i, r=beta
%                B=A+mu_1K_1C
syms s 
TT=vpa(diff(r,s),3)
T0=vpa(subs(TT,s,s_val),6)
m0=vpa(inner_s(T0,T0),6)
Tu=vpa(diff(TT,s),3)
K1=vpa(onorm(Tu),3)
K10=vpa(subs(onorm(Tu),s,s_val),3)
N11=1/K1*Tu;
N10m1=subs(N11,s,s_val)
m1=vpa(inner_s(N10m1,N10m1),6)
N1=vpa(m1*1/K1*Tu,3);
N10=vpa(subs(N1,s,s_val),6)
Tconj=co(TT);
tt=vpa(sop(N1,transpose(Tconj)),3)
tt0=vpa(subs(tt,s,s_val),6)

t=[tt(2),tt(3),tt(4),tt(5),tt(6),tt(7),tt(8)]
t0=vpa(subs(t,s,s_val),6)
e0=vpa(subs(inner(t,t),s,s_val),6)

tu=diff(t,s);
ka=vpa(sqrt(abs(inner(tu,tu))),3);
ka10=vpa(subs(ka,s,s_val),3)
%%%%calculating n1 and e1  at s_val%%%%%
n11=vpa(1/(ka)*tu,3);
e1=vpa(subs(inner(n11,n11),s,s_val),3)
n1=vpa(e1*n11,3)
n10=vpa(subs(n1,s,s_val),3)

%%%%calculating n5 and e5 at s_val %%%%%
n5=vpa(transpose(vectorproduct(t,n1)),3);
n50=vpa(subs(n5,s,s_val),3)
e5=vpa(subs(inner(n5,n5),s,s_val),3)

%%%%calculating n2 and kappa2 at s_val %%%%%
v1=vpa((diff(n1,s)-e0*(inner(diff(n1,s),t))*t-e5*(inner(diff(n1,s),n5)*n5)),3)
ka2=vpa(onorm_p(v1),3)
ka20=vpa(subs(ka2,s,s_val),3)
n22=vpa((1/ka2)*v1,6);
e2=vpa(subs(inner(n22,n22),s,s_val),6) 
n2=vpa(e2*n22,3);
n20=vpa(subs(n2,s,s_val),3)

%%%%calculating n3 at s_val %%%%%
n3=vpa(transpose(vectorproduct(n1,n2)),3)
n30=vpa(transpose(vectorproduct(n10,n20)),6)

%%%%calculating n4 at s_val %%%%%
n4=vpa(transpose(vectorproduct(n2,t)),3)
n40=vpa(transpose(vectorproduct(n20,t0)),3)

%%%%calculating n6 at s_val %%%%%
n6=vpa(transpose(vectorproduct(n2,n5)),3)
n60=vpa(transpose(vectorproduct(n20,n50)),3)

%%%%% e3-e6
e3=vpa(inner(n30,n30),6)
e4=vpa(inner(n40,n40),6)
e6=vpa(inner(n60,n60),6)
%%%%% alpha_1-alpha_6
aa10=subs(diff(n1,s),s,s_val)
alpha_1=vpa(inner(aa10,n50),6)
aa20=subs(diff(n2,s),s,s_val)
alpha_2=vpa(inner(aa20,n30),6)
alpha_3=vpa(inner(aa20,n40),6)
alpha_4=vpa(inner(aa20,n60),6)
aa30=vpa(subs(diff(n3,s),s,s_val),3)
alpha_5=vpa(inner(aa30,n40),6)
alpha_6=vpa(inner(aa30,n60),6)



tt0=vpa([0,t0(1), t0(2),t0(3), t0(4), t0(5),t0(6), t0(7)],3);
nn10=vpa([0,n10(1), n10(2),n10(3), n10(4), n10(5),n10(6), n10(7)],3);
nn20=vpa([0,n20(1), n20(2), n20(3), n20(4), n20(5),n20(6),n20(7)],3);
nn30=vpa([0,n30(1), n30(2), n30(3), n30(4), n30(5),n30(6),n30(7)],3);
nn40=vpa([0,n40(1), n40(2), n40(3), n40(4), n40(5),n40(6),n40(7)],3);
nn50=vpa([0,n50(1), n50(2),n50(3), n50(4), n50(5),n50(6),n50(7)],3);
nn60=vpa([0,n60(1), n60(2),n60(3), n60(4), n60(5),n60(6),n60(7)],3);

%%%%calculating N2 at s_val %%%%%
N20=sop(nn10,T0)

%%%%calculating N3 at s_val %%%%%
N30=sop(nn20,T0)

%%%%calculating N4 at s_val %%%%%
N40=sop(nn30,T0)

%%%%calculating N5 at s_val %%%%%
N50=sop(nn40,T0)

%%%%calculating N6 at s_val %%%%%
N60=sop(nn50,T0)

%%%%calculating N6 at s_val %%%%%
N70=sop(nn60,T0)



%%%%% m2-m7
m2=inner_s(N20,N20)
m3=inner_s(N30,N30)
m4=inner_s(N40,N40)
m5=inner_s(N50,N50)
m6=inner_s(N60,N60)
m7=inner_s(N70,N70)

Y10=m0*sop(tt0,N10)
Y20=sop(nn10,N10)
Y30=sop(nn20,N10)
Y40=sop(nn30,N10)
Y50=sop(nn40,N10)
Y60=sop(nn50,N10)
Y70=sop(nn60,N10)

y(1,1)= m0*inner_s(Y10,T0)

y(4,2)= m3*inner_s(Y20,N30)
y(5,2)= m4*inner_s(Y20,N40)
y(6,2)= m5*inner_s(Y20,N50)
y(7,2)= m6*inner_s(Y20,N60)
y(8,2)= m7*inner_s(Y20,N70)


y(5,3)= m4*inner_s(Y30,N40)
y(6,3)= m5*inner_s(Y30,N50)
y(7,3)= m6*inner_s(Y30,N60)
y(8,3)= m7*inner_s(Y30,N70)

y(6,4)= m5*inner_s(Y40,N50)
y(7,4)= m6*inner_s(Y40,N60)
y(8,4)= m7*inner_s(Y40,N70)

y(7,5)= m6*inner_s(Y50,N60)
y(8,5)= m7*inner_s(Y50,N70)

y(8,6)= m7*inner_s(Y60,N70)


m=[m0,m1,m2,m3,m4,m5,m6,m7]
C=vpa([0,0,0,0,0,0,0,0;
    y(1,1),0,0,0,0,0,0,0;
    0,0,0,y(4,2),y(5,2),y(6,2),y(7,2),y(8,2);
    0,0,-m(3)*m(4)*y(4,2),0,y(5,3),y(6,3),y(7,3),y(8,3);
    0,0,-m(3)*m(5)*y(5,2),-m(4)*m(5)*y(5,3), 0,y(6,4),y(7,4),y(8,4);
    0,0,-m(3)*m(6)*y(6,2),-m(4)*m(6)*y(6,3), -m(5)*m(6)*y(6,4),0,y(7,5),y(8,5);
    0,0,-m(3)*m(7)*y(7,2),-m(4)*m(7)*y(7,3), -m(5)*m(7)*y(7,4),-m(6)*m(7)*y(7,5),0,y(8,6);
    0,0,-m(3)*m(8)*y(8,2),-m(4)*m(8)*y(8,3), -m(5)*m(8)*y(8,4),-m(6)*m(8)*y(8,5),-m(7)*m(8)*y(8,6),0],6)

A=vpa([0,(m1*K10),0,0,0,0,0,0;
    0,0,(m0*e1*ka10),0,0,0,0,0;
    0,-(e0*ka10),0,(e2*ka20),0,0,(e5*alpha_1),0;
    0,0,(-e1*ka20),0,e3*alpha_2,e4*alpha_3,0,e6*alpha_4;
    0,0,0,-e2*alpha_2,0,e4*alpha_5,0,e6*alpha_6;
    0,0,0,-e2*alpha_3,-e3*alpha_5,0,e0*e5*ka20,-e0*e6*alpha_2;
    0,0,-e1*alpha_1,0,0,-e0*e4*ka20,0,0;
    0,0,0,-e2*alpha_4,-e3*alpha_6,e0*e4*alpha_2,0,0],6)

B=vpa(A+m1*K10*C,6)

G=vpa([T0(1),T0(2),T0(3),T0(4),T0(5),T0(6),T0(7),T0(8);
    N10(1),N10(2),N10(3),N10(4),N10(5),N10(6),N10(7),N10(8);
    N20(1),N20(2),N20(3),N20(4),N20(5),N20(6),N20(7),N20(8);
    N30(1),N30(2),N30(3),N30(4),N30(5),N30(6),N30(7),N30(8); 
    N40(1),N40(2),N40(3),N40(4),N40(5),N40(6),N40(7),N40(8);
    N50(1),N50(2),N50(3),N50(4),N50(5),N50(6),N50(7),N50(8);
    N60(1),N60(2),N60(3),N60(4),N60(5),N60(6),N60(7),N60(8);
    N70(1),N70(2),N70(3),N70(4),N70(5),N70(6),N70(7),N70(8)],6)

dGds=vpa(B*G,6) 
end

%Conjugate of the split octonion
function p=co(a)
p=vpa([a(1),-a(2),-a(3),-a(4),-a(5),-a(6),-a(7),-a(8)],6)
end

%inner product in R^7_4
function s=inner(a,b)
s=vpa(-a(1)*b(1)-a(2)*b(2)-a(3)*b(3)+a(4)*b(4)+a(5)*b(5)+a(6)*b(6)-a(7)*b(7),6);
end

%inner product in R^8_4
function s=inner_s(a,b)
s=vpa(a(1)*b(1)-a(2)*b(2)-a(3)*b(3)-a(4)*b(4)+a(5)*b(5)+a(6)*b(6)+a(7)*b(7)-a(8)*b(8),6);
end
%Norm in R^8_4
function s=onorm(a)
s=vpa(sqrt(abs(inner_s(a,a))),6)
end
%Norm in R^7_4
function s=onorm_p(a)
s=vpa(sqrt(abs(inner(a,a))),6)
end
%Split Octonion Product in R^8_4
function p=sop(a,b)
sa=a(1);
sb=b(1);
va=[a(2) a(3) a(4) a(5) a(6) a(7) a(8)];
vb=[b(2) b(3) b(4) b(5) b(6) b(7) b(8)];
v1=sa*vb;
v2=sb*va;
v3=transpose(vectorproduct(va,vb));
vi=inner(va,vb);
v=vpa(v1+v2+v3,6);
p=vpa([sa*sb-vi v(1) v(2) v(3) v(4) v(5) v(6) v(7)],6);
end
%Vector product in R^7_4
function s=vectorproduct(a,b)
a1=a(1);
a2=a(2);
a3=a(3);
a4=a(4);
a5=a(5);
a6=a(6);
a7=a(7);
b1=b(1);
b2=b(2);
b3=b(3);
b4=b(4);
b5=b(5);
b6=b(6);
b7=b(7);
A=[0 a6 -a5 -a7 a3 -a2 a4;
   -a6 0 a4 -a3 -a7 a1 a5;
   a5 -a4 0 a2 -a1 -a7 a6;
   -a7 -a3 a2 0 -a6 a5 a1;
   a3 -a7 -a1 a6 0 -a4 a2;
   -a2 a1 -a7 -a5 a4 0 a3;
   -a4 -a5 -a6 a1 a2 a3 0];
B=[b1;
b2;
b3;
b4;
b5;
b6;
b7];
s=vpa(A*B,4)
end