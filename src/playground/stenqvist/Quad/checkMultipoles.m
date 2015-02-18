function [] = checkMultipoles()
close all
clear all
clc

files_box = { 'hist_quad_xx_box.dat' 'hist_quad_xy_box.dat' 'hist_quad_xz_box.dat' 'hist_quad_yx_box.dat' 'hist_quad_yy_box.dat' 'hist_quad_yz_box.dat' 'hist_quad_zx_box.dat' 'hist_quad_zy_box.dat' 'hist_quad_zz_box.dat' };
files = { 'hist_quad_xx.dat' 'hist_quad_xy.dat' 'hist_quad_xz.dat' 'hist_quad_yx.dat' 'hist_quad_yy.dat' 'hist_quad_yz.dat' 'hist_quad_zx.dat' 'hist_quad_zy.dat' 'hist_quad_zz.dat' };

leg = {'xx' 'xy' 'xz' 'yx' 'yy' 'yz' 'zx' 'zy' 'zz'};

colors = {'k' 'r' 'b' 'r--' 'k' 'm' 'b--' 'm--' 'k' };

hFig = figure(1);
subplot(1,2,1)
maxV = zeros(1,length(files_box));
hold on
for k = 1:length(files_box)
    A = load(files_box{k});
    plot(A(:,1),A(:,2),colors{k})
    maxV(k) = max(A(:,2));
end
plot([0 0],[0 1.2*max(maxV)],'k')
legend(leg)
title('Box')
hold off

subplot(1,2,2)
maxV = zeros(1,length(files_box));
hold on
for k = 1:length(files)
    A = load(files{k});
    plot(A(:,1),A(:,2),colors{k})
    maxV(k) = max(A(:,2));
end
plot([0 0],[0 1.2*max(maxV)],'k')
legend(leg)
title('Sphere')
hold off

set(hFig, 'Position', [80 80 1400 800])



disp('Done!')