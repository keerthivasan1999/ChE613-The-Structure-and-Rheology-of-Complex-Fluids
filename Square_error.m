%ChE613A (The Structure and Rheology of Complex fluids)
%Assignment 05 (due date: 27/09/2021)
%Student name: Keerthi Vasan M (Roll No.: 21102023}


N=7;%N is number of modes in the model
%Note: If N is changed then do the corresponding change in the 'fun' statement that is used for the 'lsqcurvefit' curve fitting


%------------------------------------------
%Given datas
f=[0.0031 0.0046 0.0068 0.0099 0.0145 0.0216 0.0314 0.0462 0.0681 0.1003 0.1478 0.2122 0.3127 0.4607 0.6744 1.0000  1.4544 2.1568];
G_prime=[1.0834 1.4177 1.7969 2.2774 2.8134 3.4534 4.0792 4.8804 5.6189 6.3868 7.1673 7.9407 8.6857 9.5005 14.8752 11.5872 13.1708 15.9610];
f_doublePgiven=[0.0021 0.0031 0.0046 0.0067 0.0099 0.0146 0.0213 0.0316 0.0459 0.0677 0.0990 0.1469 0.2150 0.3167 0.46067 0.6787 0.9936 1.4544];
G_doublePgiven=[1.1330 1.3469 1.6012 1.8435 2.0955 2.3515 2.5887 2.7955 2.9236 3.0382 3.0774 3.0382 2.9995 2.9424 2.8497 2.7423 2.6389 2.6389];
%------------------------------------------



%------------------------------------------
Tau=[];
omega=f*2*3.14;
Tau_max=power(omega(1),-1);
Tau_min=power(omega(length(omega)),-1);

omega_doublePgiven=f_doublePgiven*2*3.14;
for i=1:N
    Tau(i)=Tau_min*power((Tau_max/Tau_min),((i-1)/(N-1))); 
end
%{
fun=[];
for m=1:length(omega)
    fun(m)=power((omega(m)*Tau(10)),2)/(1+power((omega(m)*Tau(10)),2));
end
plot(fun,G_prime,'o');
%}
abscissa=[];
for n=1:N
    for m=1:length(omega)
        abscissa(n,m)=power((omega(m)*Tau(n)),2)/(1+power((omega(m)*Tau(n)),2));  
    end
end
x0=[];
for n=1:N
    x0(n)=1;
end
options = optimoptions('lsqcurvefit','StepTolerance',exp(pi/2));
lb = [];
ub = [];
fun=@(x,abscissa)(x(1)*abscissa(1,:))+(x(2)*abscissa(2,:))+(x(3)*abscissa(3,:))+(x(4)*abscissa(4,:))++(x(5)*abscissa(5,:))+(x(6)*abscissa(6,:))+(x(7)*abscissa(7,:));
RM=lsqcurvefit(fun,x0,abscissa,G_prime,lb,ub,options);
for i=1:N
    fprintf('Mode no.: %i, G_i = %f, Tau_i = %f \n',i,RM(i),Tau(i))
end
G_dis="G=";
Eta_dis="Eta=";
for i=1:N
    if i==N
        G_dis=G_dis+"{"+num2str(RM(i))+"*exp(-t/"+num2str(Tau(i))+")}";
        Eta_dis=Eta_dis+"{"+num2str(RM(i)*Tau(i))+"*exp(-t/"+num2str(Tau(i))+")}";
    else 
        G_dis=G_dis+"{"+num2str(RM(i))+"*exp(-t/"+num2str(Tau(i))+")} + ";
        Eta_dis=Eta_dis+"{"+num2str(RM(i)*Tau(i))+"*exp(-t/"+num2str(Tau(i))+")} + ";
    end
end
fprintf('Relaxation modulus is \n%s',G_dis)
fprintf('\nShear  viscosity is \n%s',Eta_dis)
%{

x=[];
x0=[];
for n=1:N
    x0(n)=1;
    fun_sub=@(x,abscissa)(x(n)*abscissa(n,:));
end
fun=fun_sub;
x=lsqcurvefit(fun,x0,abscissa,G_prime);
%}
%G_doubleprime=[];
sum=0;
for i=1:length(omega)
    for j=1:N
        sum=sum+((RM(j)*omega(i)*Tau(j))/(1+power((omega(i)*Tau(j)),2)));
    end
    G_doubleprime(i)=sum;
    sum=0;
end

%disp(G_doubleprime);
%plot(omega,G_doubleprime,'o',omega_doublePgiven,G__doublePgiven,'.');
%loglog(f,G_prime,'o',f,G_doubleprime,'.')
error=[];
for i=1:length(omega)
    error(i)=power((G_doubleprime(i)-G_doublePgiven(i)),2)/G_doublePgiven(i);
end
disp(error)
% G=0;
% t=1:1000;
% for i=1:N
%     G=G+(RM(i)*exp(-t/Tau(i)) );
% end
% plot(t,G);
% title('Plot of relaxation modulus G(t) with respect to time')
% xlabel('Time (s)')
% ylabel('G(t)')
