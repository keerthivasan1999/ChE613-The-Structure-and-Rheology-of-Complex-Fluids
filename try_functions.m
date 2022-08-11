RM=[3.3144 5.6098 2.5231 2.5965 1.5966 1.6556];
Tau=[0.0738 0.2734 1.0121 3.7471 13.8736 51.3663];
t=1:100;
G=0;
for i=1:length(RM)
    G=G+(RM(i)*exp(-t/Tau(i)) );
end
%plot(t,G);
G_dis="G=";
for i=1:length(RM)
    if i==length(RM)
        G_dis=G_dis+"{"+num2str(RM(i))+"*exp(-t/"+num2str(Tau(i))+")}";
    else 
        G_dis=G_dis+"{"+num2str(RM(i))+"*exp(-t/"+num2str(Tau(i))+")} + ";
    end
end
fprintf('Relaxation modulus is \n%s',G_dis)