
function [xv,yv] = get_xv_and_yv(image)
    [x,y] = find(image);
    x = x / 100 - 1; y = y/100 -1;
    n = size(x);
    
    angle = [];
    for i=1:n
        angle(end+1) = atan2(y(i),x(i));
    end
    %sort x and y in decreasing angles:
    [angle, idx] = sort(angle, 'descend');
    x = x(idx);
    x(end+1) = x(1);
    y = y(idx); 
    y(end+1) = y(1);
      
    xv = (x + 1);
    yv = (y + 1);
end