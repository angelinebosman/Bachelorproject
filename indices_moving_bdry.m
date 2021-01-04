function [placement, edge] = indices_moving_bdry(i,j,deltax, r1, xc, yc)
    x = deltax*(i-1/2);
    y = deltax*(j-1/2);
    n= size(xc,1);
    T = linspace(0,2*pi,n);
    xv = r1 + r1*cos(T); xc = reshape(xc, [1,n]);
    xv = [xv NaN xc];
    yv = r1 + r1*sin(T);  yc = reshape(yc, [1,n]);
    yv = [yv NaN yc];
    if inpolygon(x,y,xv,yv) == 1
        x_plus = deltax*(i+1/2);x_min = deltax*(i-3/2);
        y_plus = deltax*(j+1/2);y_min = deltax*(j-3/2);
        if inpolygon(x_plus,y,xv,yv) == 0
            if x > r1 %punt zit aan de rechterkant van het domein
                placement = "buitenrand";
                if inpolygon(x,y_min,xv,yv) == 0
                    edge = "west_noord"; %
                elseif inpolygon(x,y_plus,xv,yv) == 0
                    edge = "oost_noord";
                else
                    edge = "noord"; %
                end
            else %punt ligt aan de linkerkant van het domein
                placement = "binnenrand";
                if inpolygon(x,y_min,xv,yv) == 0
                    edge =  "oost_zuid";
                elseif inpolygon(x,y_plus,xv,yv) == 0
                    edge =  "west_zuid";
                else
                    edge = "zuid";
                end
            end
        elseif inpolygon(x_min,y,xv,yv) == 0
            if x > r1
                placement = "binnenrand";
                if inpolygon(x,y_min,xv,yv) == 0
                    edge = "oost_noord";
                elseif inpolygon(x,y_plus,xv,yv) == 0
                    edge = "west_noord";
                else
                    edge = "noord"; 
                end
            else
                placement = "buitenrand";
                if inpolygon(x,y_min,xv,yv) == 0
                    edge = "west_zuid";
                elseif inpolygon(x,y_plus,xv,yv) == 0
                    edge = "oost_zuid";
                else
                    edge = "zuid"; %
                end
            end
        elseif inpolygon(x,y_min,xv,yv) == 0
            if y > r1
                placement = "binnenrand";
                if inpolygon(x_min,y,xv,yv) == 0
                    edge = "oost_zuid"
                elseif inpolygon(x_plus,y,xv,yv) == 0
                    edge = "oost_noord"
                else
                    edge = "oost"; 
                end
            else
                placement = "buitenrand";
                if inpolygon(x_min,y,xv,yv) == 0
                    edge = "west_zuid";
                elseif inpolygon(x_plus,y,xv,yv) == 0
                    edge = "west_noord";
                else
                    edge = "west"; %
                end
            end
        elseif inpolygon(x,y_plus, xv, yv) == 0
            if y > r1
                placement = "buitenrand";
                if inpolygon(x_min,y,xv,yv) == 0
                    edge = "oost_zuid";
                elseif inpolygon(x_plus,y,xv,yv) == 0
                    edge = "oost_noord"; %
                else
                    edge = "oost"; %
                end
            else
                placement = "binnenrand";
                if inpolygon(x_min,y,xv,yv) == 0
                    edge = "zuid_west";
                elseif inpolygon(x_plus,y,xv,yv) == 0
                    edge = "zuid_oost";
                else
                    edge = "west";
                end
            end
        else
            placement = "inside";
            edge = NaN;
        end
    else
        placement = "outside";
        edge = NaN;
    end
end