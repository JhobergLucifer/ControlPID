%# Control PID rootlocus
%# by Jhoberg Quevedo Ruiz
%# Engineer


%           u    +         +---+      +------+y
%          ------>(+)----->| Gs |----->| Gp |--------+------->
%                  ^ -     +---+      +------+       |
%                  |                                 |
%                  +---------------------------------+


pkg load symbolic;
pkg load control;
pkg load signal;

syms S
syms Kp
syms Kd
syms Ki
syms E
syms wn
syms t

Kpow=49.5
s=tf('s')
  
Gpt=Kpow*exp(-t/14.940)+19
Gps=Kpow*exp(-t/14.940)+19
Gp=laplace(Gpt)
Gp=simplify(Gp)
GP=Kpow*-((S-(2/0.3))/(S+(2/0.3)))/(0.53*S+1)
Gptf=(102339*s + 1900)/(2*s*(747*s + 50))
GPTF=(102339*S + 1900)/(2*S*(747*S + 50))

figure(1)
hold off
t1 = linspace(0, 100, 1000);  % vector time
y1=Kpow*exp(-t1/14.940)+19
plot(t1,y1)
title ("Step Temperature");
xlabel ("Time");
ylabel ("Temperature");
%text (pi, 0.7, "arbitrary text");
legend ("Signal temperature To");



%timeTs=0.1
%E=5%
%Td=0%
  
[numGp, denGp] = numden(sym(Gp))
numGp=expand(numGp)
denGp=expand(denGp)
den=double(coeffs(denGp))
num=double(coeffs(numGp))

Gp1=tf(Gptf)
[z,p,k] = tf2zp(num,den)
[p,z]=pzmap(Gp1)
  
  %set figure 
h = findobj(gca, 'type', 'line');
set(h, 'markersize', 9)
text(real(roots(num)) - 0.1, imag(roots(num)) + 0.1, 'Zero')
text(real(roots(den)) - 0.1, imag(roots(den)) + 0.1, 'Pole')
axis equal

figure(2)

step(Gp1);
title ("Step Temperature");
xlabel ("Time");
ylabel ("Temperature");
%text (pi, 0.7, "arbitrary text");
legend ("Signal temperature To");

%step(feedback(Gp1))
%title ("Step Temperature");
%xlabel ("Time");
%ylabel ("Temperature");
%text (pi, 0.n7, "arbitrary text");
legend ("Signal Temperature To");
hold off



%%%%%%%%%%%%%%%%Control PID

%error state stable rampa unitaria <=0,001
%sobrepaso max<=0.01%
%tiempo levantamiento tr<=0,001 seg
%tiempo asentamiento ts <=0,007 seg




%calculo PID controller
%1) error estado estable ess
%Ke=limS->0 Gol(S), entrada escalon unitario
 limit(GP, S, 0)
  Ke=49.5
  
%ess=1/(1+Ke)
solve(0.001==1/(1+Ki),Ki)
%Kp=100
  
%%%%%%%%
%2)caculo Overshot sistema 2nd order
clear E
syms E
solve(0.001==exp(-pi*E/(1-E*E)^0.5),E)

  E=0.9102;

%%%%%
%3) calculo Tiempo retardo tr 
% aproxlinea recta , fr=1/(L*C)^0.5=316.22 , tr>>10*1/fr
clear wn
syms wn
tr=(0.8+2.5*E)/(wn)
wn=solve(0.001==(0.8+2.5*E)/wn,wn)
wn=3075.4

%4) calculo Tiempo ascentamiento ts
  Ts=4/(wn*E)
  Ts=1.4286e-03
       
%5)caculo PID , ecuacion caracteristica Close Loop Gcl
%divisor orden 2 , S*S+2*E+wn+wn*wn

sd1=-E*wn+i*wn*(1-E*E)^0.5
sd2=-E*wn-i*wn*(1-E*E)^0.5

sd1=-2799.2 + 1273.9i
     
     
%Gc=K*(s+a)*(s+b)/s
     

%sum(Z)-sum(P)=-180
    % (teta1-teta4)-(teta2-teta3)=-180
     teta4=180-atan(1273.9/(2799.2-0))
     teta3=180-atan(1273.9/(2799.2-0.01856))
     teta3=180-atan(1273.9/(2799.2-0.066934))
     %teta1=atan(182.73/(a-134.68))
%     solve(atan(42/(a-40))-(179.26+179.18)=-180,a)
     %teta1=atan(1766/a-3075.7)
%sum(tetaZ)-sum(tetaP)=-180
     a=1273.9/tan(-0.43)+2799.2
     a=21.525

%Koverall=IILP/IILZ
%Kall=L2*L3/L1
%     L2=(42*42+38*38)^0.5
%     L3=(42*42+39*39)^0.5
%     L1=(42*42+(43.66-40)*(43.66-40))^0.5
%    Kall=L2*L3/L1
     
% |Gc(S)*H(S)|=1     
% K*(S+43.66)*1/((S+1)*(S+2))=1
%
%K=((S+1)*(S+2))/(S+43.66)
     KHs=(1/GPTF)*(1/(S+21.525))
S=abs(sd1)
K=eval(KHs)

K = 0.014498
     


%ess
 limit(GP, S, 0)
Ke=49.5
%Kc=K/Ke
 Kc=K/Ke
     
%PID final
%Gc=Kp+Ki/S+Kd*S
%Gc=(Kp*S+Ki+Kd*S²)/S
%Gc=Kd*(S²+Kp/Kd*S+Ki/Kd)/S
     
%Gc=K*(S+a)*(S+b)/S
%Gc=K*(S²+S*a+S*b+a*b)/S
%Gc=Kd*(S*S+87.32*S+4.366)/S    

clear S
syms S     
Gc=0.014498*(S+21.525)*(S+0.1)/S
Gc=expand(Gc)
     
Kd1=0.01691
Kp1=0.3657
Ki1=0.036399

%Kd*X=0.0638
     %X=0.3657/0.016910
Kp=21.626
Ki=2.1525
Kd=0.014498
     
Gc1=Kp+'s'*Kd+Ki/'s'
%Gc2=tf([Kd Kp Ki],[1 0])


 
Gc1=tf([Kd Kp Ki],[1 0])
Gcl1=Gc1*Gp1/(1+Gc1*Gp1)
figure(3)
     hold on
step(Gcl1,0.1)

     figure(3)
     hold off
title ("PID Control Temperature");
xlabel ("Time");
ylabel ("Temperature");
%text (pi, 0.7, "arbitrary text");
legend ("Signal temperature To");


   
%ajuste erros stado stable ess
%rampa S*GP
limit(S*GP, S, 0)
ts=4/(wn*E)
ts=1.9820E-3

figure(4)
subplot (2, 1, 1)
step(Gp1);
title ("Step Temperature Model Physical");
xlabel ("Time");
ylabel ("Temperature");
%text (pi, 0.7, "arbitrary text");
legend ("Signal temperature To");

%step(feedback(Gp1))
%title ("Step Temperature Model Physical");
%xlabel ("Time");
%ylabel ("Temperature");
%text (pi, 0.n7, "arbitrary text");
legend ("Signal Temperature To");


subplot (2, 1, 2)
step(Gcl1,0.1)
title ("Step Temperature PID Control ");
xlabel ("Time");
ylabel ("Temperature");
%text (pi, 0.7, "arbitrary text");
legend ("Signal temperature To");

%step(feedback(Gp1))
%title ("Step Temperature PID Controller");
%xlabel ("Time");
%ylabel ("Temperature");
%text (pi, 0.n7, "arbitrary text");
legend ("Signal Temperature To");

hold off

