load('data_LD_2.mat', 'data');
x = reshape(data(:, 1), [128, 1001]);
t = reshape(data(:, 2), [128, 1001]);
n = reshape(data(:, 3), [128, 1001]);
u = reshape(data(:, 4), [128, 1001]);
p = reshape(data(:, 5), [128, 1001]);
q = reshape(data(:, 6), [128, 1001]);
E = reshape(data(:, 7), [128, 1001]);

colormap('hot')
subplot(5, 1, 1);
pcolor(t, x, n);
ylabel('x(\pi/k_1)');
shading interp;  % 使用插值方法填充颜色
cb = colorbar;
cb.Label.String = 'n'; % 这里设置颜色条的标签名

subplot(5, 1, 2);
pcolor(t, x, u);
ylabel('x(\pi/k_1)');
shading interp;  % 使用插值方法填充颜色
cb = colorbar;
cb.Label.String = 'u'; % 这里设置颜色条的标签名

subplot(5, 1, 3);
pcolor(t, x, p);
ylabel('x(\pi/k_1)');
shading interp;  % 使用插值方法填充颜色
cb = colorbar;
cb.Label.String = 'p'; % 这里设置颜色条的标签名

subplot(5, 1, 4);
pcolor(t, x, q);
ylabel('x(\pi/k_1)');
shading interp;  % 使用插值方法填充颜色
cb = colorbar;
cb.Label.String = 'q'; % 这里设置颜色条的标签名

subplot(5, 1, 5);
pcolor(t, x, E);
ylabel('x(\pi/k_1)');
xlabel('t(\omega_{pe}^{-1})');
shading interp;  % 使用插值方法填充颜色
cb = colorbar;
cb.Label.String = 'E'; % 这里设置颜色条的标签名


width = 640; % 画布宽度
height = 480; % 画布高度
set(gcf, 'Position', [100, 100, width, height]);
