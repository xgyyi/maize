% ���ú�����
  % KS��Kennard-Stone�㷨����������
  % LDA��LDA�㷨�����ݽ�ά
  % KNN: K�����㷨����LDA���ݽ���Ʒ���б�ģ��
  
% �������ݴ���
  % Writing by Guiyan yang
  % 2021.11.2

clear all
close all
clc
%% �����������ӹ�������
fea = xlsread('ABCDE_Ave.xlsx');% �������ӹ���
gnd = [ones(35,1);2*ones(35,1);3*ones(35,1);4*ones(35,1);5*ones(35,1)];% ����ǩ

% ����Kennard-Stone�㷨������ѵ�����Ͳ��Լ���������
N = 116;% ѵ������������2:1����
[fea_Train,fea_Test,trainIdx,dminmax,testIdx] = KS(fea,N);
gnd_Train = gnd(trainIdx);
gnd_Test = gnd(testIdx);

% ����Ϊѵ������������λ��
Train_1 = find(gnd_Train==1);
Train_2 = find(gnd_Train==2);
Train_3 = find(gnd_Train==3);
Train_4 = find(gnd_Train==4);
Train_5 = find(gnd_Train==5);
%% ����ѵ����LDAͶӰ������ͶӰֵ
% ��������
options = [];
options.Fisherface = 1;
% ����LDA�㷨
[eigvector_Train, eigvalue_Train] = LDA(gnd_Train, options, fea_Train);
% ����LDAͶӰֵ
LDA_Train = fea_Train*eigvector_Train;  
%% ��������ƽ����������ͼ
x = xlsread('X.xlsx');% ������
In = 1:20:451;
figure(1)
h1 = plot(x,mean(fea(1:35,:))+160,'k-o','MarkerIndices', In);
hold on
h2 = plot(x,mean(fea(36:70,:))+120,'k-*','MarkerIndices', In);
h3 = plot(x,mean(fea(71:105,:))+80,'k-^','MarkerIndices', In);
h4 = plot(x,mean(fea(106:140,:))+40,'k-d','MarkerIndices', In);
h5 = plot(x,mean(fea(141:175,:))-7,'k-s','MarkerIndices', In);
xlim([380 1850])
ylim([-0.05 430])
legend([h1(1) h2(1) h3(1) h4(1) h5(1)],'֣��958','ũ��86','ũ��87','�Ǻ�605','����968','fontsize',8,'FontName','SimSun')
xlabel('���� / cm^{-1} ','fontsize',8,'FontName','SimSun')
ylabel('ǿ��','fontsize',8,'FontName','SimSun','Rotation',90)
set(gca,'FontName','Times New Roman','FontSize',8)
set(gcf,'unit','centimeters','position',[10 5 15 8]) 
%% ����ѵ�������������б����ͼ
% ֣��958
figure(21)
h1 = plot(LDA_Train(Train_1,1),LDA_Train(Train_1,2),'ko');
legend(h1(1),'֣��958','fontsize',8,'FontName','SimSun')
xlabel('��һ������� ','fontsize',8,'FontName','SimSun')
ylabel('�ڶ��������','fontsize',8,'FontName','SimSun','Rotation',90)
set(gca,'FontName','Times New Roman','FontSize',8)
set(gcf,'unit','centimeters','position',[10 5 12 8]) 

% 5����������
figure(22)
h1 = plot(LDA_Train(Train_1,1),LDA_Train(Train_1,2),'ko');
hold on
h2 = plot(LDA_Train(Train_2,1),LDA_Train(Train_2,2),'k*');
h3 = plot(LDA_Train(Train_3,1),LDA_Train(Train_3,2),'k^');
h4 = plot(LDA_Train(Train_4,1),LDA_Train(Train_4,2),'kd');
h5 = plot(LDA_Train(Train_5,1),LDA_Train(Train_5,2),'ks');
legend([h1(1) h2(1) h3(1) h4(1) h5(1)],'֣��958','ũ��86','ũ��87','�Ǻ�605','����968','fontsize',8,'FontName','SimSun','Box','off')
xlabel('��һ������� ','fontsize',8,'FontName','SimSun')
ylabel('�ڶ��������','fontsize',8,'FontName','SimSun','Rotation',90)
set(gca,'FontName','Times New Roman','FontSize',8)
set(gcf,'unit','centimeters','position',[10 5 12 8]) 
%% ����Ʒ���б�ģ��
% ������Լ�LDAͶӰֵ
LDA_Test = fea_Test*eigvector_Train;  
% ����K�����㷨������ģ��
[rate,Result]= KNN(LDA_Train,gnd_Train,LDA_Test,gnd_Test,1);% ���Լ�rate ׼ȷ�ʣ�Result label ���









