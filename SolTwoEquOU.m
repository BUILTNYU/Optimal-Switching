function [ equ_diff_1, equ_diff_2] = SolTwoEquOU(Q_h,Q_l)
%input:
% Q_h=15;
% Q_l=10;

l=0.2; %Mean-reverting speed parameter
sigma=7; % Constant of volatility
rou=0.07; % interest rate
m=40; % mean value of demand 
% miu=0.0; %dfit
a1=rou/2/l; %kummer function A parameter1
a2=0.5; % kunmmer function A parameter 2
% a3=l/(sigma^2)*(m-Q)^2; % kummer function A parameter 3
b1=0.5*(1+rou/l); % kummer function B parameter 1
b2=1.5; % kummer function B parameter 2
% b3=l/sigma^2*(m-Q)^2;  % kummer function B parameter 3
f_plus=10;
f_minus=10;
w1=6.378;%para
% w2=185.51;

delta0=-l^0.5/sigma*gamma(0.5)*gamma(0.5*(1+rou/l))/gamma(1.5)/gamma(rou/2/l); % relationship between A0 and B0
delta1=l^0.5/sigma*gamma(0.5)*gamma(0.5*(1+rou/l))/gamma(1.5)/gamma(rou/2/l); % relationship between A0 and B0

delta_capQ=0.001;
% gE_deri_capQh=(get_E_capQ(Q_h+delta_capQ)-get_E_capQ(-delta_capQ+Q_h))/2/delta_capQ;
% gE_deri_capQl=(get_E_capQ(Q_l+delta_capQ)-get_E_capQ(-delta_capQ+Q_l))/2/delta_capQ;
gE_deri_capQh=w1/(rou+l);
gE_deri_capQl=w1/(rou+l);

gHah_deri_capQ=(kummerCal(a1,a2,l/(sigma^2)*(m-(Q_h+delta_capQ))^2)-kummerCal(a1,a2,l/(sigma^2)*(m-(Q_h-delta_capQ))^2))/2/delta_capQ;

gHal_deri_capQ=(kummerCal(a1,a2,l/(sigma^2)*(m-(Q_l+delta_capQ))^2)-kummerCal(a1,a2,l/(sigma^2)*(m-(Q_l-delta_capQ))^2))/2/delta_capQ;

gHbh_deri_capQ=(kummerCal(b1,b2,l/(sigma^2)*(m-(Q_h+delta_capQ))^2)-kummerCal(b1,b2,l/(sigma^2)*(m-(Q_h-delta_capQ))^2))/2/delta_capQ;

gHbl_deri_capQ=(kummerCal(b1,b2,l/(sigma^2)*(m-(Q_l+delta_capQ))^2)-kummerCal(b1,b2,l/(sigma^2)*(m-(Q_l-delta_capQ))^2))/2/delta_capQ;

G0h=gHah_deri_capQ+delta0*(m-Q_h)*gHbh_deri_capQ-delta0*kummerCal(b1,b2,l/(sigma^2)*(m-Q_h)^2);
G1h=gHah_deri_capQ+delta1*(m-Q_h)*gHbh_deri_capQ-delta1*kummerCal(b1,b2,l/(sigma^2)*(m-Q_h)^2);
G0l=gHal_deri_capQ+delta0*(m-Q_l)*gHbl_deri_capQ-delta0*kummerCal(b1,b2,l/(sigma^2)*(m-Q_l)^2);
G1l=gHal_deri_capQ+delta1*(m-Q_l)*gHbl_deri_capQ-delta1*kummerCal(b1,b2,l/(sigma^2)*(m-Q_l)^2);

Hah=kummerCal(a1,a2,l/(sigma^2)*(m-Q_h)^2);
Hal=kummerCal(a1,a2,l/(sigma^2)*(m-Q_l)^2);

Hbh=kummerCal(b1,b2,l/(sigma^2)*(m-Q_h)^2);
Hbl=kummerCal(b1,b2,l/(sigma^2)*(m-Q_l)^2);

K0h=Hah+delta0*(m-Q_h)*Hbh;
K1h=Hah+delta1*(m-Q_h)*Hbh;

K0l=Hal+delta0*(m-Q_l)*Hbl;
K1l=Hal+delta1*(m-Q_l)*Hbl;

equ_diff_1=((get_E_capQ(Q_h)-f_plus)*G0h-K0h*gE_deri_capQh)/(K0h*G1h-K1h*G0h)-(K0l*gE_deri_capQl-(get_E_capQ(Q_l)+f_minus)*G0l)/(K1l*G0l-K0l*G1l);
equ_diff_2=((get_E_capQ(Q_h)-f_plus)*G1h*K0h-K0h*K1h*gE_deri_capQl)/K0h/(K0h*G1h-K1h*G0h)-(K1l*K0l*gE_deri_capQl-(get_E_capQ(Q_l)+f_minus)*G1l*K0l)/K0l/(K1l*G0l-G1l*K0l);

end