global it

F = getframe(gcf);
Im = frame2im(F);
[A, map] = rgb2ind(Im, 256);
if it == 10
    % 删除旧的 GIF 文件（如果存在）
    if exist('phase.gif', 'file')
        delete('phase.gif');
    end
    % 第一次写入时创建 GIF 文件，并确保使用 GIF89a 格式
    imwrite(A, map, 'phase.gif', 'gif', 'Loopcount', Inf, 'DelayTime', 0.1);
else
    % 后续追加帧
    imwrite(A, map, 'phase.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
end
