% 参数定义
h0 = 485; % 工作行程 
t=8; % 放大系数
theta = 0.1; % 修正角系数
process_rate = 4; % 生产率
delta0 = pi; % 推程角  
cut1 = pi / 3; % 远休角
cut2 = pi / 3; % 近休角

% 推程回程压力角
alpha1 = 40 * pi / 180;  
alpha2 = 70 * pi / 180;

omega = (process_rate * 360 / 60) * pi / 180; % 将角速度转换为弧度每秒

delta1 = theta * delta0;
delta2 = theta * delta0;

%行程
h = h0 / t;
h1 = (delta1 * h)/(2*delta0-delta1-delta2);
h2 = (delta2 * h)/(2*delta0-delta1-delta2);

delta3 = 2 * pi - cut1 - cut2 - delta0; % 回程角

% 计算步长
delta_all = linspace(0, 2 * pi, 3600);

% 分割区间
s = zeros(size(delta_all));
v = zeros(size(delta_all));
a = zeros(size(delta_all));
ds = zeros(size(delta_all));
r0 = zeros(size(delta_all));

for i = 1:length(delta_all)
    delta = delta_all(i);
    if delta <= delta1
        s(i) = h1 * ((delta / delta1) - sin(pi * delta / delta1) / pi);
        v(i) = h1 * omega * (1 - cos(pi * delta / delta1)) / delta1;
        a(i) = pi * h1 * omega^2 * sin(pi * delta / delta1) / delta1^2;
        ds(i) = h1 / delta1 * (1 - cos(pi * delta / delta1));
        r0(i) = abs(ds(i) / tan(alpha1) - s(i));
    elseif delta <= delta0 - delta2
        s(i) = h1 + (h - h1 - h2) * (delta - delta1) / (delta0 - delta1 - delta2);
        v(i) = ((h - h1 - h2) * omega) / (delta0 - delta1 - delta2);
        a(i) = 0;
        ds(i) = (h - h1 - h2) / (delta0 - delta1 - delta2);
        r0(i) = abs(ds(i) / tan(alpha1) - s(i));
    elseif delta <= delta0
        s(i) = h - h2 * ((delta0 - delta) / delta2) + h2 * sin(pi * (delta0 - delta) / delta2) / pi;
        v(i) = h2 * omega / delta2 - h2 * omega * cos(pi * (delta0 - delta) / delta2) / delta2;
        a(i) = -h2 * omega^2 * pi * sin(pi * (delta0 - delta) / delta2) / delta2^2;
        ds(i) = h2 / delta2 * (1 - cos(pi * (delta0 - delta) / delta2));
        r0(i) = abs(ds(i) / tan(alpha1) - s(i));
    elseif delta <= delta0 + cut1
        s(i) = h;
        v(i) = 0;
        a(i) = 0;
        ds(i) = 0;
        r0(i) = NaN;  % 休止没有r0
    elseif delta <= delta0 + cut1 + delta3
        s(i) = h * (1 - (delta - delta0 - cut1) / delta3 + (1/(2*pi)) * sin(2*pi * (delta - delta0 - cut1) / delta3));
        v(i) = (h * omega / delta3) * (cos(2 * pi * (delta - delta0 - cut1) / delta3) - 1);
        a(i) = -2 * h * pi * omega^2 / delta3 * sin(2 * pi * (delta - delta0 - cut1) / delta3);
        ds(i) = -h / (delta3) + h * (1/(2*pi)) * cos((2*pi * (delta - delta0 - cut1)) / delta3) * (2*pi / (delta3));
        r0(i) = abs(abs(ds(i) / tan(alpha2)) - s(i));
    else
        s(i) = 0;
        v(i) = 0;
        a(i) = 0;
        ds(i) = 0;
        r0(i) = NaN;  % 休止没有r0
    end
end

% 输入最大基圆半径
max_r0 = max(r0(~isnan(r0)));
prompt = {'计算所得的最大基圆半径为 ' + string(max_r0) + '，你希望使用什么数值计算？'};
dlgtitle = '输入新的基圆半径';
dims = [1 50];
definput = {'60.625'};  % 设置默认值为60.625
answer = inputdlg(prompt, dlgtitle, dims, definput);
R_base = str2double(answer{1});  % 更新基圆半径

% 更新理论和实际轮廓线半径
R_theoretical = R_base + s;
R_actual = R_theoretical - 0.3 * R_base;  % 实际轮廓线半径减去滚子半径

% 创建表格并写入Excel
T = table(delta_all' * 180 / pi, s', v', a', ds', r0', R_theoretical', R_actual', ...
          'VariableNames', {'角度 (Degrees)', '位移 (s)', '速度 (v)', '加速度 (a)', '导数 (ds/delta)', '基圆半径 r_0', '理论廓线', '实际廓线'});
filename = '水平凸轮.xlsx';
writetable(T, filename, 'Sheet', 1);
disp(['Data written to ', filename]);

% 设定颜色
color1 = [0.00, 0.45, 0.74];  % 蓝色
color2 = [0.85, 0.33, 0.10];  % 红色
color3 = [0.47, 0.67, 0.19];  % 绿色
color4 = [0.64, 0.08, 0.18];  % 深红色

% 图像生成和布局
figure;
tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% 绘制位移曲线
nexttile
plot(delta_all * 180/pi, s, '-', 'LineWidth', 2, 'Color', color1);
title('位移曲线 s(\delta)');
xlabel('\delta (度)');
ylabel('位移 s (mm)');
grid on;
legend('位移', 'Location', 'best');

% 绘制速度曲线
nexttile
plot(delta_all * 180/pi, v, '-', 'LineWidth', 2, 'Color', color2);
title('速度曲线 v(\delta)');
xlabel('\delta (度)');
ylabel('速度 v (mm/s)');
grid on;
legend('速度', 'Location', 'best');

% 绘制加速度图像
nexttile
plot(delta_all * 180/pi, a, '-', 'LineWidth', 2, 'Color', color3);
title('加速度曲线 a(\delta)');
xlabel('\delta (度)');
ylabel('加速度 a (mm/s^2)');
grid on;
legend('加速度', 'Location', 'best');

% 绘制导数函数图像
nexttile
plot(delta_all * 180/pi, ds, '-', 'LineWidth', 2, 'Color', color4);
title('导数函数图像 ds/d\delta');
xlabel('\delta (度)');
ylabel('导数 ds/d\delta (mm/度)');
grid on;
legend('导数', 'Location', 'best');

% 绘制基圆半径 r0 曲线
nexttile
plot(delta_all * 180/pi, r0, '-', 'LineWidth', 2, 'Color', color1);
title('基圆半径曲线 r_0(\delta)');
xlabel('\delta (度)');
ylabel('基圆半径 r_0 (mm)');
grid on;
legend('基圆半径', 'Location', 'best');

% 绘制极坐标图表示理论和实际轮廓线
nexttile
polarplot(delta_all, R_theoretical, '-', 'LineWidth', 2, 'Color', color1);
hold on;
polarplot(delta_all, R_actual, '-', 'LineWidth', 2, 'Color', color2);
title('凸轮轮廓');
legend('理论轮廓线', '实际轮廓线');
grid on;

% 调整图形窗口的大小
set(gcf, 'Position', [100, 100, 1200, 800]);