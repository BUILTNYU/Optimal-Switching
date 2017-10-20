x=20;
X=[x];
delta=1;
sigma=2.1;
driftRate=0.01;
for i=0:240
x=x+driftRate*delta+sigma*sqrt(delta)*normrnd(0,0.6);
if x<0|x>100
    break;
end
X=[X;x];
end
y=0:1:length(X)-1;
y=y*delta
plot(y,X)
xlabel('Month')
ylabel('Demand density (trips/mile^2/hr)')