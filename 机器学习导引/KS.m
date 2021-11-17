% Kennard-Stone算法
% 把所有的样本都看作训练集候选样本，依次从中挑选样本进训练集：
% 首先选择欧氏距离最远的两个向量对进入训练集；
% 在接下来的迭代过程中，分别计算剩余的样本与已选择的两个样本之间的距离，拥有最大最小距离的待选样本被选入训练库．以此类推，达到所要求的样本数目。
% 该方法优点：能保证训练库中样本按照空间距离分布均匀。
% 缺点：需要进行数据转换和计算样本两两空间距离，计算量大。

function [X_Selected,X_Rest,m,dminmax,NotSelectedSample] = KS(X,N)

% Kennard-Stone Algorithm for selection of samples
% [m,dminmax] = ks(X,N);
%
% X --> Matrix of instrumental responses，每行为一个样本，列为变量
% N --> Number of samples to be selected (minimum of 2)，为期望的训练集样本数
%
% X_Selected --> 挑选出的训练集数据阵，其与m中的序号对应
% X_Rest --> 为挑选出的测试集数据阵
% m --> Indexes of the selected samples，为挑选的训练集样本编号，行向量，1*N
% NotSelectedSample -->选择的预测集样本序号，升序排列，行向量，1*N
%
% dminmax(1) = 0;第一个值为0
% dminmax(2) = Euclidean distance between the two first samples selected by the algorithm
%              欧式距离，通过算法选择的前两个样本
% dminmax(i) = Smallest distance between the i-th selected sample and the previously selected ones (i > 2)
%              最小的距离，第i个选择的样本与先前选择的之间的最小距离

%%
%%********************* start of the kennard-stone step one **********************%%
%　1.首先计算两两样本之间的距离
dminmax = zeros(1,N); % Initializes the vector of minimum distances 
M = size(X,1); % Number of rows in X (samples)
samples = 1:M;

D = zeros(M,M); % Initializes the matrix of distances 距离矩阵
for i=1:M-1
    xa = X(i,:);
    for j = i+1:M
      xb = X(j,:);
      D(i,j) = norm(xa - xb);% 计算的是欧氏距离，norm(A)或者norm(A,2)计算的就是A的2范数；2范数：
      % 计算步骤是先计算A*A'（这里A'代表转置，也就是原矩阵*（原矩阵的转置）），然后计算他们乘积的特征值，取最大的那个特征值开根号即可
    end
end

% D: Upper Triangular Matrix 上三角矩阵
% D(i,j) = Euclidean distance between objects i and j (j > i)

%　选择距离最大的两个样品
[maxD,index_row] = max(D); % 对每一列求最大值，得到一行最大值向量，最大值所在的行， maxD = Row vector containing the largest element of each column in D
                             % index_row(n) = Index of the row with the largest element in the n-th column

[dummy,index_column] = max(maxD); % index_column = column corresponding to the largest element in matrix D

m(1) = index_row(index_column);
m(2) = index_column;

dminmax(2) = D(m(1),m(2));
%%************************ end of the kennard-stone step one *****************************%%

%%
%　2.然后分别计算剩余的样本与已选择的两个样本之间的距离,对于每个剩余样本而言, 其与已
%　选样品之间的最短距离被选择,然后选择这些最短距离中相对最长的距离所对应的样本, 作
%　为第三个样品,重复步骤,直至所选的样品的个数等于事先确定的数目为止%%********************* start of the kennard-stone step two **********************%%

for i=3:N
    % This routine determines the distances between each sample still available for selection and each of the samples already selected
    pool = setdiff(samples,m); % c = setdiff(A,B) 返回在A中有，而B中没有的值，结果向量将以升序排序返回。pool = Samples still available for selection 
    dmin = zeros(1,M-i+1); % Initializes the vector of minimum distances between each sample in pool and the samples already selected
    for j = 1:(M-i+1) % For each sample xa still available for selection
        indexa = pool(j); % indexa = index of the j-th sample in pool (still available for selection)

        d = zeros(1,i-1); % Initializes the vector of distances between the j-th sample in pool and the samples already selected
        for k = 1:(i-1) % The distance with respect to each sample already selected is analyzed
            indexb =  m(k); % indexb = index of the k-th sample already selected
            if indexa < indexb
                d(k) = D(indexa,indexb);
            else
                d(k) = D(indexb,indexa);
            end
        end
        dmin(j) = min(d);
    end
    % The selected sample corresponds to the largest dmin
    [dminmax(i),index] = max(dmin);
    m(i) = pool(index);
end
%%************************ end of the kennard-stone step two *****************************%%

%%
%%************************ start of export the result *****************************%%
%　反馈输出结果
SelectedRowIndex = m; 

for i=1:length(m) 
    
    X_Selected(i,:) = X(SelectedRowIndex(i),:); 
end 

NotSelectedSample = setdiff(samples,m); 
for i=1:length(NotSelectedSample) 
    
    X_Rest(i,:) = X(NotSelectedSample(i),:); 
end 

%%************************ end of export the result *******************************%%















