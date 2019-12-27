clear all
clc
num=[0.5];
den=[2 4 5 1];
H=tf(num,den);
H=feedback(H,1);
tempo=[0:0.05:100];
[amostra]=step(H,tempo);
a=[tempo' amostra];
csvwrite('csvlist.csv',a);

hold on
plot(a(:,1),a(:,2),'g')
step(H)

