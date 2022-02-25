#############################
#                           #
# Tarefa Computacional 1    #
# Bruno Hayashi - 214012    #
#                           #
#############################
clear all; clc; close all;
hold on;
##### PROBLEMA 1 #####

%TAREFA 1.1
    ##DADOS
    L=3;                      %[m]
    E=200e9;                  %[N/m^2]
    A=20e-3.^2;               %[m^2]
    q0=10e3;                  %[N/m]
    q= q0;                     %[N/m]
    ua=0;                     %[m]
    ub=[0,0.001,0.01];        %[m]
    
    p=0.1; %passo
    x=[0:p:L];
    ##Solucao Analitica
    ux11= @(i,x) ((q*L)/(E*A*2)*(L*x-x.^2))+ub(i)*x/L;
    
    for i=1:3
      n=0;
      for xn=0:p:3
        n=n+1;
        S11(i,n)=ux11(i,xn);
      endfor
    endfor
    ##Carregamento
    n=0;
    for xn=0:p:3;
      n=n+1;
      q11(n)=q;
    endfor
    ##Graficos
    subplot(2,1,2),plot(x,S11(1,:),x,S11(2,:),x,S11(3,:));
    legend(['ub= ' num2str(ub(1))],['ub= ' num2str(ub(2))],['ub= ' num2str(ub(3))]);
    subplot(2,1,1),plot(x,q11);
%TAREFA 1.2
    ##DADOS
    L=3;                      %[m]
    E=200e9;                  %[N/m^2]
    A=20e-3.^2;               %[m^2]
    q0=10e3;                  %[N/m]
    q= @(x) q0*x/L;           %[N/m]
    ua=0;                     %[m]
    ub=[0,0.001,0.01];        %[m]
    
    p=0.1; %passo
    x=[0:p:L];
    ##Solucao Analitica
    ux12= @(i,x) q0/(6*E*A)*(x*L-(x^3/L))+ub(i)*x/L;
    
    for i=1:3
      n=0;
      for xn=0:p:3
        n=n+1;
        S12(i,n)=ux12(i,xn);
      endfor
    endfor
    ##Carregamento
    n=0;
    for xn=0:p:3;
      n=n+1;
      q12(n)=q(xn);
    endfor
    ##Graficos
    figure(2)
    subplot(2,1,2),plot(x,S12(1,:),x,S12(2,:),x,S12(3,:));
    legend(['ub= ' num2str(ub(1))],['ub= ' num2str(ub(2))],['ub= ' num2str(ub(3))]);
    subplot(2,1,1),plot(x,q12);
%TAREFA 1.3
    ##DADOS
    L=3;                      %[m]
    E=200e9;                  %[N/m^2]
    A=20e-3.^2;               %[m^2]
    q0=10e3;                  %[N/m]
    q= @(x) q0*sin(pi*x/L);           %[N/m]
    ua=0;                     %[m]
    ub=[0,0.001,0.01];        %[m]
    
    p=0.01; %passo
    x=[0:p:L];
    ##Solucao Analitica
    ux13= @(i,x) L.^2*q0/(pi.^2*E*A)*sin(pi*x/L)+ub(i)*x/L;
    for i=1:3
      n=0;
      for xn=0:p:3
        n=n+1;
        S13(i,n)=ux13(i,xn);
      endfor
    endfor
    ##Carregamento
    n=0;
    for xn=0:p:3;
      n=n+1;
      q13(n)=q(xn);
    endfor
    ##Graficos
    figure(3)
    subplot(2,1,2),plot(x,S13(1,:),x,S13(2,:),x,S13(3,:));
    legend(['ub= ' num2str(ub(1))],['ub= ' num2str(ub(2))],['ub= ' num2str(ub(3))]);
    subplot(2,1,1),plot(x,q13);

%TAREFA 1.4
    ##DADOS
    L=3;                      %[m]
    E=200e9;                  %[N/m^2]
    A=20e-3.^2;               %[m^2]
    q0=10e3;                  %[N/m]
    q= @(x) q0*x/L;           %[N/m]
    ua=0;                     %[m]
    ub=0;                     %[m]
    n=40;                     %Numero de Pontos
    p=L/(n-1);          %passo
    ##VETORES E MATRIZES
    xi=[0:p:L];
    xi(1)=[]; xi(end)=[]; %Correcao do vetor x
    k= -q0*p^2/(L*E*A)*xi;
    k(1)=k(1)-ua; k(end)=k(end)-ub; %Ajuste
    M=zeros(n-2,n-2);
    for i=2:n-2
      M(i,i-1)=1;
      M(i,i)=-2;
      M(i,i+1)=1;
    endfor
    M(1,1)=-2;M(1,2)=1; M(:,end)=[]; %Ajuste
    
    ##Resolvendo equacao matricial
    u=M\k';
    u=[ua;u;ub];
    x=[0,xi,L];
    ##Solucao Analitica
    uA= q0/(6*E*A)*(x*L-(x.^3/L))+ub*x/L;
    TA=q0/(6*A)*(L-3*x.^2/L)+ub*E/L;
    ##Calculo das tensoes
    for i=2:n-1
      T(i)=E*(u(i+1)-u(i-1))/(2*p);
    endfor
    T(1)=E*(u(2)-u(1))/p;
    T(n)=E*(u(n)-u(n-1))/p;
    ##Erro
    erro=((uA-u'));
    normaErro=norm(erro);
    normaRelativa=(norm(erro)/norm(uA))*100;
    erroT=((TA-T));
    normaErroT=norm(erroT);
    normaRelativaT=(norm(erroT)/norm(TA))*100;
    ##Graficos
    figure(4)
    subplot(2,2,1),plot(x,u,'-ro',x,uA,'g'); grid on;
    legend('numérico','analítico');
    subplot(2,2,3),plot(x,erro); grid on;
    legend([ 'Norm.Rel % =' num2str(normaRelativa)])
    title(['Norma do Erro N = ' num2str(n)]);
    xlabel(' Posição x');
    ylabel(' Erro=(analítico-numérico) ');
    subplot(2,2,2),plot(x,T,x,TA)
    legend('numérico','analítico')
    subplot(2,2,4),plot(x,erroT);
%TAREFA 1.5
    ##DADOS
    L=3;                      %[m]
    E=200e9;                  %[N/m^2]
    A=20e-3.^2;               %[m^2]
    q0=10e3;                  %[N/m]
    q= @(x) q0*sin(pi*x/L);   %[N/m]
    ua=0;                     %[m]
    ub=0;                     %[m]
    n=100;                     %Numero de Pontos
    p=L/(n-1);          %passo
    ##VETORES E MATRIZES
    xi=[0:p:L];
    xi(1)=[]; xi(end)=[]; %Correcao do vetor x
    k= -q0*p^2*sin(pi*xi/L)/(E*A);
    k(1)=k(1)-ua; k(end)=k(end)-ub; %Ajuste
    M=zeros(n-2,n-2);
    for i=2:n-2
      M(i,i-1)=1;
      M(i,i)=-2;
      M(i,i+1)=1;
    endfor
    M(1,1)=-2;M(1,2)=1; M(:,end)=[]; %Ajuste
    
    ##Resolvendo equacao matricial
    u=M\k';
    u=[ua;u;ub];
    x=[0,xi,L];
    ##Solucao Analitica
    uA= L.^2*q0/(pi.^2*E*A)*sin(pi*x/L)+ub*x/L;
    TA=L*q0/(pi*E*A)*cos(pi*x/L)+ub/L;
    ##Calculo das tensoes
    for i=2:n-1
      T(i)=E*(u(i+1)-u(i-1))/(2*p);
    endfor
    T(1)=[];
    ##Erro
    erro=abs((uA-u'));
    normaErro=norm(erro);
    normaRelativa=(norm(erro)/norm(uA))*100;
    ##Graficos
    figure(5)
    subplot(2,2,1),plot(x,u,'-ro',x,uA,'g'); grid on;
    legend('numérico','analítico');
    subplot(2,2,3),plot(x,erro); grid on;
    legend([ 'Norm.Rel % =' num2str(normaRelativa)])
    title(['Norma do Erro N = ' num2str(n)]);
    xlabel(' Posição x');
    ylabel(' Erro=(analítico-numérico) ');
    subplot(2,2,2),plot(xi,T,x,TA)
    legend('numérico','analítico')
%TAREFA 2.1
    ##DADOS
    L=10;                      %[m]
    E=210e9;                  %[N/m^2]
    A=20e-3.^2;               %[m^2]
    q0=10e3;                  %[N/m]
    q= @(x) q0*sin(pi*x/L);   %[N/m]
    ua=0;                     %[m]
    ub=0;                     %[m]
    n=5;                     %Numero de Pontos
    p=L/(n-1);          %passo
    
    hl=40e-3;
    h0=3*hl;
    F=100e3;
    b=30e-3;
    ##VETORES E MATRIZES
    x=[0:p:L];

    
    ##Solucao analitica
    uA=F*L/(2*hl*E*b)*(log(3*hl)-log(3*hl-2*hl*x/L));
    
    

    ## solucao numerica
    ## E*(A(x)u"(x)+A'(x)u'(x))=0
    A= b*(3*hl-2*hl*x/L); %função da área
    for i=2:n-1
      dA(i)=(A(i+1)-A(i-1))/(2*p); %primeira derivada da área (MDF centrada)
    endfor
    dA(1)=(A(2)-A(1))/p;           %área em x=0 (MDF adiantada)
    dA(n)=(A(n)-A(n-1))/p;         %área em x=l (MDF atrasada)
    k=zeros(1,n+1);                %vetor da direita
    k(1)=k(1)-ua;                  %condição de contorno de Dirichlet em x=0
    k(end)=F*p/(E*A(n));           %Condição de contorno de Neuman em x=L
                                   %Constante definidas descrevendo u' e u" com MDF centrada
    K1=@(i) -dA(i)*p^2+A(i)*2*p;   %Constante 1 da matriz
    K2=@(i)-A(i)*4*p;              %Constante 2
    K3=@(i)dA(i)*p^2+A(i)*2*p;     %Constante 3
    M=zeros(n,n);
    
    for i=2:n                      %Colocando as constante na diagonal da matriz
      M(i,i-1)=K1(i);
      M(i,i)=K2(i);
      M(i,i+1)=K3(i);
    endfor
    M(1,1)=K2(1);M(1,2)=K3(1);     %Ajuste da primeira linha da matriz
    M(end+1,:)=zeros(1,n+1);       %Adicionando uma linha nula no final para deixar quadrada
    M(end,end-1)=-1; M(end,end)=1; %Deixando a última linha descrevendo condição de contorno (MDF atrasado)
    
    ##resolvendo
    u=M\k';
    u(end)=[];                     %Último valor não tem significado (u(n+1))
    
    plot(x,uA,x,u)
    legend('analítico','numérico')
    ##Erro
    erro=abs((uA-u'));
    normaErro=norm(erro);
    normaRelativa=(norm(erro)/norm(uA))*100;