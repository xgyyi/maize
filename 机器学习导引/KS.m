% Kennard-Stone�㷨
% �����е�����������ѵ������ѡ���������δ�����ѡ������ѵ������
% ����ѡ��ŷ�Ͼ�����Զ�����������Խ���ѵ������
% �ڽ������ĵ��������У��ֱ����ʣ�����������ѡ�����������֮��ľ��룬ӵ�������С����Ĵ�ѡ������ѡ��ѵ���⣮�Դ����ƣ��ﵽ��Ҫ���������Ŀ��
% �÷����ŵ㣺�ܱ�֤ѵ�������������տռ����ֲ����ȡ�
% ȱ�㣺��Ҫ��������ת���ͼ������������ռ���룬��������

function [X_Selected,X_Rest,m,dminmax,NotSelectedSample] = KS(X,N)

% Kennard-Stone Algorithm for selection of samples
% [m,dminmax] = ks(X,N);
%
% X --> Matrix of instrumental responses��ÿ��Ϊһ����������Ϊ����
% N --> Number of samples to be selected (minimum of 2)��Ϊ������ѵ����������
%
% X_Selected --> ��ѡ����ѵ��������������m�е���Ŷ�Ӧ
% X_Rest --> Ϊ��ѡ���Ĳ��Լ�������
% m --> Indexes of the selected samples��Ϊ��ѡ��ѵ����������ţ���������1*N
% NotSelectedSample -->ѡ���Ԥ�⼯������ţ��������У���������1*N
%
% dminmax(1) = 0;��һ��ֵΪ0
% dminmax(2) = Euclidean distance between the two first samples selected by the algorithm
%              ŷʽ���룬ͨ���㷨ѡ���ǰ��������
% dminmax(i) = Smallest distance between the i-th selected sample and the previously selected ones (i > 2)
%              ��С�ľ��룬��i��ѡ�����������ǰѡ���֮�����С����

%%
%%********************* start of the kennard-stone step one **********************%%
%��1.���ȼ�����������֮��ľ���
dminmax = zeros(1,N); % Initializes the vector of minimum distances 
M = size(X,1); % Number of rows in X (samples)
samples = 1:M;

D = zeros(M,M); % Initializes the matrix of distances �������
for i=1:M-1
    xa = X(i,:);
    for j = i+1:M
      xb = X(j,:);
      D(i,j) = norm(xa - xb);% �������ŷ�Ͼ��룬norm(A)����norm(A,2)����ľ���A��2������2������
      % ���㲽�����ȼ���A*A'������A'����ת�ã�Ҳ����ԭ����*��ԭ�����ת�ã�����Ȼ��������ǳ˻�������ֵ��ȡ�����Ǹ�����ֵ�����ż���
    end
end

% D: Upper Triangular Matrix �����Ǿ���
% D(i,j) = Euclidean distance between objects i and j (j > i)

%��ѡ���������������Ʒ
[maxD,index_row] = max(D); % ��ÿһ�������ֵ���õ�һ�����ֵ���������ֵ���ڵ��У� maxD = Row vector containing the largest element of each column in D
                             % index_row(n) = Index of the row with the largest element in the n-th column

[dummy,index_column] = max(maxD); % index_column = column corresponding to the largest element in matrix D

m(1) = index_row(index_column);
m(2) = index_column;

dminmax(2) = D(m(1),m(2));
%%************************ end of the kennard-stone step one *****************************%%

%%
%��2.Ȼ��ֱ����ʣ�����������ѡ�����������֮��ľ���,����ÿ��ʣ����������, ������
%��ѡ��Ʒ֮�����̾��뱻ѡ��,Ȼ��ѡ����Щ��̾����������ľ�������Ӧ������, ��
%��Ϊ��������Ʒ,�ظ�����,ֱ����ѡ����Ʒ�ĸ�����������ȷ������ĿΪֹ%%********************* start of the kennard-stone step two **********************%%

for i=3:N
    % This routine determines the distances between each sample still available for selection and each of the samples already selected
    pool = setdiff(samples,m); % c = setdiff(A,B) ������A���У���B��û�е�ֵ��������������������򷵻ء�pool = Samples still available for selection 
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
%������������
SelectedRowIndex = m; 

for i=1:length(m) 
    
    X_Selected(i,:) = X(SelectedRowIndex(i),:); 
end 

NotSelectedSample = setdiff(samples,m); 
for i=1:length(NotSelectedSample) 
    
    X_Rest(i,:) = X(NotSelectedSample(i),:); 
end 

%%************************ end of export the result *******************************%%















