% 奇异谱分析的Matlab实现
function filteredSignal = SSA(signal,L,R)
    %% 参数
    %————输入————
    % signal为一维时间序列
    % L为滑动窗口长度，L∈[2,N/2]
    % R为奇异值数量(即取前R个最大的奇异值进行趋势项提取)，R∈[1,L]
    %————输出————
    % filteredSignal为SSA提取得到的有用信号
    N=length(signal);
    if(L>N/2)
        L=N-L;
    end
    if(R>L)
        R=L;
    end
    %% 组建轨迹矩阵
    K=N-L+1;
    X=zeros(L,K);
    for i=1:L
        for j=1:K
            X(i,j)=signal(j+i-1);
        end
    end
    %% 奇异值分解SVD
    Xt=X';
    S=X*Xt;
    [U,lambda]=eig(S);%eig返回矩阵的特征值和特征向量，U是特征向量，lambda是特征值
    [~,order]=sort(diag(lambda),'descend');
    U=U(:,order);%根据特征值大小对特征向量进行排序
    V=Xt*U;
    %% 分组（只取了前R个较大的奇异值对应的信号分量进行还原）
    Y=U(:,1:R)*V(:,1:R)';
    %% 对角平均
    filteredSignal=zeros(N,1);
    Kp=max(L,K);
    Lp=min(L,K);
    for k=1:Lp
        for p=1:k
            filteredSignal(k,1)=filteredSignal(k,1)+Y(p,k-p+1)/k;
        end
    end
    for k=Lp+1:Kp
        for p=1:Lp
            filteredSignal(k,1)=filteredSignal(k,1)+Y(p,k-p+1)/Lp;
        end
    end
    for k=Kp+1:N
        for p=k-Kp+1:N-Kp+1
            filteredSignal(k,1)=filteredSignal(k,1)+Y(p,k-p+1)/(N-k+1);
        end
    end
end