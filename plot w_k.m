%%%%% to plot the w-k with the data from solving the cold dispersion relationship
%%%%%Jingchun LI, E-mail:lijingchun2016@gmail.com
clc;clear;
z=load(['data.dat']);
KK=z(1,:);
x=sort(z(2:11,:),1,'descend');
WW=x(1:5,:);

 h1=plot(KK,WW(4,:),'LineWidth',2);hold on;
 h2=plot(KK,WW(5,:),'LineWidth',2);hold on;
axis([0,20,0,30])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mp=1.6726231e-27;me=mp/10^10;epsilon0=8.854187817E-12;
qe=1.60217662e-19;q_i=qe;q_e=-qe;B0=1;m_e=me;m_i=mp;
n_e=225*epsilon0*B0^2/mp;n_i=n_e;
wci=q_i*B0/m_i;
wpi=sqrt(n_i*q_i^2/(epsilon0*m_i));
wce=abs(q_e*B0/m_e);
wlh=1/sqrt(1/(wci^2+wpi^2)+1/(wce*wci));
wlh=wlh/wci;
wlh=wlh*ones(280);
%wce=wce/wci;
%wce=wce*ones(280);
h3=plot(KK,wlh,'g--','LineWidth',2);
%h2=plot(KK,wce,'y--','LineWidth',2);
legend([h3],'\Omega_{lh}','Location','SouthEast');
ylabel('\omega/\Omega_p');xlabel('k\lambda_p');set(gca,'FontSize',15);

