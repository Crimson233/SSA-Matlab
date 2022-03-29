load("Signal.mat")

L = 1800;%滑动窗口大小
R=4;%提取的奇异值数
filteredSignal=SSA(Signal,L,R);

figure
plot(Signal,'k');
hold on
plot(filteredSignal,'r','LineWidth',2);
legend("原序列",'提取的趋势分量');