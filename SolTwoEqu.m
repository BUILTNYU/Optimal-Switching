function [ equ_diff_1, equ_diff_2] = SolTwoEqu( Q_h,Q_l)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%variables to find
% Q_h=15;
% Q_l=10;


% l=1; %Mean-reverting speed parameter
sigma=9; % Constant of volatility
rou=0.06; % interest rate
% m=40; % mean value of demand 
miu=0.05; % grow rate

f_plus=10;
f_minus=10;
w1=6.378;%para
%meter in Et(Expectation)
w2=185.51;
%Et=(1/(rou-miu)*(w1*Q_h)+w2/rou);

c1=0.5*((1-2*miu/sigma^2)+((1-2*miu/sigma^2)^2+8*rou/sigma^2)^0.5);
c2=0.5*((1-2*miu/sigma^2)-((1-2*miu/sigma^2)^2+8*rou/sigma^2)^0.5);

P1=((1-c1)*w1/(rou-miu)+c1/Q_h*(f_plus+w2/rou))/Q_h^(c2-1);
P2=(w1*(1-c1)/(rou-miu)+c1/Q_l*(w2/rou-f_minus))/Q_l^(c2-1);
equ_diff_1=P1-P2;

P3=(c2*(w2/rou+f_plus)+Q_h*w1*(1-c2)/(rou-miu))/Q_h^c1;
P4=(c2*(w2/rou-f_minus)+Q_l*(1-c2)*w1/(rou-miu))/Q_l^c1;
equ_diff_2=P3-P4;
end




