% 调用函数：
  % KS：Kennard-Stone算法，样本划分
  % LDA：LDA算法，数据降维
  % KNN: K近邻算法，对LDA数据建立品种判别模型
  
% 样本数据处理
  % Writing by Guiyan yang
  % 2021.11.2

clear all
close all
clc
%% 导入玉米种子光谱数据
fea = xlsread('ABCDE_Ave.xlsx');% 玉米种子光谱
gnd = [ones(35,1);2*ones(35,1);3*ones(35,1);4*ones(35,1);5*ones(35,1)];% 类别标签

% 调用Kennard-Stone算法，划分训练集和测试集样本划分
N = 116;% 训练集样本数，2:1划分
[fea_Train,fea_Test,trainIdx,dminmax,testIdx] = KS(fea,N);
gnd_Train = gnd(trainIdx);
gnd_Test = gnd(testIdx);

% 划分为训练集的样本，位置
Train_1 = find(gnd_Train==1);
Train_2 = find(gnd_Train==2);
Train_3 = find(gnd_Train==3);
Train_4 = find(gnd_Train==4);
Train_5 = find(gnd_Train==5);
%% 计算训练集LDA投影向量及投影值
% 参数定义
options = [];
options.Fisherface = 1;
% 调用LDA算法
[eigvector_Train, eigvalue_Train] = LDA(gnd_Train, options, fea_Train);
% 计算LDA投影值
LDA_Train = fea_Train*eigvector_Train;  
%% 绘制种子平均拉曼光谱图
x = xlsread('X.xlsx');% 横坐标
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
legend([h1(1) h2(1) h3(1) h4(1) h5(1)],'郑单958','农大86','农大87','登海605','京科968','fontsize',8,'FontName','SimSun')
xlabel('波数 / cm^{-1} ','fontsize',8,'FontName','SimSun')
ylabel('强度','fontsize',8,'FontName','SimSun','Rotation',90)
set(gca,'FontName','Times New Roman','FontSize',8)
set(gcf,'unit','centimeters','position',[10 5 15 8]) 
%% 绘制训练集样本线性判别分析图
% 郑单958
figure(21)
h1 = plot(LDA_Train(Train_1,1),LDA_Train(Train_1,2),'ko');
legend(h1(1),'郑单958','fontsize',8,'FontName','SimSun')
xlabel('第一分类变量 ','fontsize',8,'FontName','SimSun')
ylabel('第二分类变量','fontsize',8,'FontName','SimSun','Rotation',90)
set(gca,'FontName','Times New Roman','FontSize',8)
set(gcf,'unit','centimeters','position',[10 5 12 8]) 

% 5种玉米种子
figure(22)
h1 = plot(LDA_Train(Train_1,1),LDA_Train(Train_1,2),'ko');
hold on
h2 = plot(LDA_Train(Train_2,1),LDA_Train(Train_2,2),'k*');
h3 = plot(LDA_Train(Train_3,1),LDA_Train(Train_3,2),'k^');
h4 = plot(LDA_Train(Train_4,1),LDA_Train(Train_4,2),'kd');
h5 = plot(LDA_Train(Train_5,1),LDA_Train(Train_5,2),'ks');
legend([h1(1) h2(1) h3(1) h4(1) h5(1)],'郑单958','农大86','农大87','登海605','京科968','fontsize',8,'FontName','SimSun','Box','off')
xlabel('第一分类变量 ','fontsize',8,'FontName','SimSun')
ylabel('第二分类变量','fontsize',8,'FontName','SimSun','Rotation',90)
set(gca,'FontName','Times New Roman','FontSize',8)
set(gcf,'unit','centimeters','position',[10 5 12 8]) 
%% 构建品种判别模型
% 计算测试集LDA投影值
LDA_Test = fea_Test*eigvector_Train;  
% 调用K近邻算法，构建模型
[rate,Result]= KNN(LDA_Train,gnd_Train,LDA_Test,gnd_Test,1);% 测试集rate 准确率；Result label 结果









