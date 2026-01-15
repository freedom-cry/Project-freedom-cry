% 参数定义
bottle_R = 60; % 瓶子大端直径
process_rate = 4; % 生产率
e = 8; % 偏置
theta = 0.1; % 修正角系数
delta0 = 180 * pi / 180;  % 推程角
cut1 = pi / 3; % 远休角
cut2 = pi / 3; % 近休角
beta_gamma = 1; % 默认为1

omega = (process_rate * 360 / 60) * pi / 180; % 将角速度转换为弧度每秒

delta1 = theta * delta0;
delta2 = theta * delta0;

delta3 = 2 * pi - cut1 - cut2 - delta0; % 回程角

alpha1 = 40 * pi / 180; % 推程回程压力角
alpha2 = 70 * pi / 180;

delta4 = cut2; % 回程角（水平凸轮远休角）
delta5 = cut1; % 推程角（水平凸轮近休角）
h = bottle_R * 1.2 / 2; 

% 定义delta的取值范围
delta = linspace(0, 2 * pi, 3600);

% 初始化数组
s = zeros(size(delta));
v = zeros(size(delta));
a = zeros(size(delta));
ds = zeros(size(delta));
r0 = zeros(size(delta));

% 计算每个区间的位移、速度、加速度、导数和 r0
for i = 1:length(delta)
    if delta(i) <= delta0
        s(i) = h;
        v(i) = 0;
        a(i) = 0;
        ds(i) = 0;
        r0(i) = NaN;
    elseif delta(i) <= delta0 + delta4
        s(i) = h * (1 - (delta(i) - delta0) / delta4 + (1/(2*pi)) * sin(2*pi * (delta(i) - delta0) / delta4));
        v(i) = (h * omega / delta4) * (cos(2 * pi * (delta(i) - delta0) / delta4) - 1);
        a(i) = -(2 * h * pi * omega^2 / delta4^2) * sin(2*pi * (delta(i) - delta0) / delta4);
        ds(i) = -h / delta4 + h * (1/(2 * pi)) * 2 * pi / delta4 * cos(2 * pi * (delta(i) - delta0) / delta4);
        r0(i) = sqrt(((abs(ds(i) - beta_gamma * e) / tan(alpha2)) - s(i))^2 + e^2);
    elseif delta(i) <= delta0 + delta4 + delta3
        s(i) = 0;
        v(i) = 0;
        a(i) = 0;
        ds(i) = 0;
        r0(i) = NaN;
    else
        s(i) = h * ((delta(i) - (delta0 + delta4 + delta3)) / delta4 - (1/(2*pi)) * sin(2*pi * (delta(i) - (delta0 + delta4 + delta3)) / delta4));
        v(i) = (h * omega / delta4) * (1 - cos(2 * pi * (delta(i) - (delta0 + delta4 + delta3)) / delta4));
        a(i) = (2 * h * pi * omega^2 / delta4) * sin(2 * pi * (delta(i) - (delta0 + delta4 + delta3)) / delta4);
        ds(i) = h / delta4 - h * (1/(2 * pi)) * 2 * pi / delta4 * cos(2 * pi * (delta(i) - (delta0 + delta4 + delta3)) / delta4);
        r0(i) = sqrt(((abs(ds(i) - beta_gamma * e) / tan(alpha1)) - s(i))^2 + e^2);
    end
end

% 计算最大基圆半径
max_r0 = max(r0(~isnan(r0)));
prompt = {'计算所得的最大基圆半径为 ' + string(max_r0) + '，你希望使用什么数值计算？'};
dlgtitle = '输入新的基圆半径';
dims = [1 50];
definput = {'60.625'};  % 设置默认值为60.625
answer = inputdlg(prompt, dlgtitle, dims, definput);

% 检查用户是否输入了值
if ~isempty(answer)
    R_base = str2double(answer{1});
else
    R_base = 60.625;  % 如果用户取消输入，则使用默认值
end

% 更新理论和实际轮廓线半径
R_theoretical = R_base + s;
R_actual = R_theoretical - 0.3 * R_base;  % 实际轮廓线半径减去滚子半径

% 创建表格并写入Excel
T = table(delta' * 180 / pi, s', v', a', ds', r0', R_theoretical', R_actual', ...
          'VariableNames', {'角度 (Degrees)', '位移 (s)', '速度 (v)', '加速度 (a)', '导数 (ds/delta)', 'r0', '理论廓线', '实际廓线'});
filename = '垂直凸轮.xlsx';
writetable(T, filename, 'Sheet', 1);
disp(['Data written to ', filename]);

%图像生成%

% 图像生成和布局
figure;
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% 绘制位移曲线
nexttile
plot(delta * 180 / pi, s, 'b-', 'LineWidth', 2);
title('位移曲线 s(\delta)');
xlabel('\delta (度)');
ylabel('位移 s (mm)');
grid on;
legend('位移', 'Location', 'best');

% 绘制速度曲线
nexttile
plot(delta * 180 / pi, v, 'r-', 'LineWidth', 2);
title('速度曲线 v(\delta)');
xlabel('\delta (度)');
ylabel('速度 v (mm/s)');
grid on;
legend('速度', 'Location', 'best');

% 绘制加速度图像
nexttile
plot(delta * 180 / pi, a, 'g-', 'LineWidth', 2);
title('加速度曲线 a(\delta)');
xlabel('\delta (度)');
ylabel('加速度 a (mm/s^2)');
grid on;
legend('加速度', 'Location', 'best');

% 绘制导数函数图像
nexttile
plot(delta * 180 / pi, ds, 'm-', 'LineWidth', 2);
title('导数函数图像 ds/d\delta');
xlabel('\delta (度)');
ylabel('导数 ds/d\delta (mm/度)');
grid on;
legend('导数', 'Location', 'best');

% 绘制基圆半径 r0 曲线
nexttile
plot(delta * 180 / pi, r0, 'k-', 'LineWidth', 2);
title('基圆半径曲线 r_0(\delta)');
xlabel('\delta (度)');
ylabel('基圆半径 r_0 (mm)');
grid on;
legend('基圆半径', 'Location', 'best');

% 绘制极坐标图表示理论和实际轮廓线
nexttile
polarplot(delta, R_theoretical, 'b-', 'LineWidth', 2, 'DisplayName', '理论轮廓线');  % 理论轮廓线
hold on;
polarplot(delta, R_actual, 'r-', 'LineWidth', 2, 'DisplayName', '实际轮廓线');  % 实际轮廓线
title('凸轮轮廓');
legend('show');
grid on;

% 调整图形窗口的大小
set(gcf, 'Position', [100, 100, 1200, 800]);

% 增强图形的视觉效果
set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'LineWidth', 1.5);