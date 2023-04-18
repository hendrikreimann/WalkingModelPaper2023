function plotArrowhead(ax, anchor, target, length, width, color)
    % transform axes coordinates to pixel coordinates
    axes_origin_pixel = getpixelposition(ax);
    xlimit = get(ax, 'xlim');
    ylimit = get(ax, 'ylim');
    
    anchor_pixel = axesToPixel(anchor);
    target_pixel = axesToPixel(target);
    v_pixel = target_pixel - anchor_pixel;
    angle = atan2(v_pixel(2), v_pixel(1));
    
    A = anchor;
    A_pixel = axesToPixel(A);
    D_pixel_x = A_pixel(1) + length * cos(angle);
    D_pixel_y = A_pixel(2) + length * sin(angle);

    B_pixel_x = A_pixel(1) + width * sin(angle);
    B_pixel_y = A_pixel(2) - width * cos(angle);
    C_pixel_x = A_pixel(1) - width * sin(angle);
    C_pixel_y = A_pixel(2) + width * cos(angle);
    
    B = pixelToAxes([B_pixel_x, B_pixel_y]);
    C = pixelToAxes([C_pixel_x, C_pixel_y]);
    D = pixelToAxes([D_pixel_x, D_pixel_y]);
    
    plot(ax, [A(1) B(1) D(1) C(1) A(1)], [A(2) B(2) D(2) C(2) A(2)])
    patch(ax, 'XData', [B(1) D(1) C(1) B(1)], 'YData', [B(2) D(2) C(2) B(2)], 'FaceColor', color, 'linewidth', 0.5, 'EdgeColor', color')

function point_pixel = axesToPixel(point_axes)
    x_pixel = axes_origin_pixel(1) + axes_origin_pixel(3) * (point_axes(1)-xlimit(1))/(xlimit(2)-xlimit(1));
    y_pixel = axes_origin_pixel(2) + axes_origin_pixel(4) * (point_axes(2)-ylimit(1))/(ylimit(2)-ylimit(1));
    point_pixel = zeros(size(point_axes));
    point_pixel(1) = x_pixel;
    point_pixel(2) = y_pixel;
end

function point_axes = pixelToAxes(point_pixel)
    w = xlimit(2) - xlimit(1);
    h = ylimit(2) - ylimit(1);
    x_axes = xlimit(1) + w * (point_pixel(1)-axes_origin_pixel(1))/axes_origin_pixel(3);
    y_axes = ylimit(1) + h * (point_pixel(2)-axes_origin_pixel(2))/axes_origin_pixel(4);
    
    point_axes = zeros(size(point_pixel));
    point_axes(1) = x_axes;
    point_axes(2) = y_axes;
end



end