function [placement, edge] = indices(i,j,deltax, r1, r2)
    x = deltax*(i-1/2);
    y = deltax*(j-1/2);
    [T,R] = meshgrid(linspace(0,2*pi,64),linspace(r1,r2,2));
    xv = r1 + R.*cos(T); xv = [xv(1,:) NaN flip(xv(2,:))];
    yv = r1 + R.*sin(T); yv = [yv(1,:) NaN flip(yv(2,:))];
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