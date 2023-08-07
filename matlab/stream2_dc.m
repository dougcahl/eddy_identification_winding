function [X,Y] = stream2_dc(x,y,u,v,startx,starty,dr,N)
% bilinear stream function
if nargin < 7
    dr = 0.1;
    N = 10000;
end
% row x
% column y
% u(x,y)
X = nan(N,1);
Y = nan(N,1);

%%
% figure
% quiver(x,y,u,v)
% hold on
% plot(startx,starty,'bs')
X(1) = startx;
Y(1) = starty;
dx = x - startx;
dy = y - starty;
d = dx.^2+dy.^2;
[~,j] = min(d(:));
[row,col] = ind2sub(size(x),j(1));

for i = 2:N
    
%     toc
%     disp('a')
%     tic
    
    x1 = x(row,col);
    y1 = y(row,col);
    u1 = u(row,col);
    v1 = v(row,col);
    if startx > x1
        row2 = row + 1;
    else
        row2 = row - 1;
    end
    if starty > y1
        col2 = col + 1;
    else
        col2 = col - 1;
    end
    x2 = x(row2,col);
    y2 = y(row2,col);
    x3 = x(row,col2);
    y3 = y(row,col2);
    x4 = x(row2,col2);
    y4 = y(row2,col2);

    u2 = u(row2,col);
    v2 = v(row2,col);
    u3 = u(row,col2);
    v3 = v(row,col2);
    u4 = u(row2,col2);
    v4 = v(row2,col2);

    r1 = sqrt( (startx-x1).^2 + (starty-y1).^2 );
    r2 = sqrt( (startx-x2).^2 + (starty-y2).^2 );
    r3 = sqrt( (startx-x3).^2 + (starty-y3).^2 );
    r4 = sqrt( (startx-x4).^2 + (starty-y4).^2 );
    if r1 == 0
        u0 = u1;
        v0 = v1;
    elseif r2 == 0
        u0 = u2;
        v0 = v2;
    elseif r3 == 0
        u0 = u3;
        v0 = v3;
    elseif r4 == 0
        u0 = u4;
        v0 = v4;
    else

        if isnan(u1) && isnan(u2) && isnan(u3) && isnan(u4)
            break
        end
        
        if isnan(u1)
            u1r1 = 0;
            v1r1 = 0;
            nm1 = 0;
        else
            u1r1 = u1/r1;
            v1r1 = v1/r1;
            nm1 = 1/r1;
        end

        if isnan(u2)
            u2r2 = 0;
            v2r2 = 0;
            nm2 = 0;
        else
            u2r2 = u2/r2;
            v2r2 = v2/r2;
            nm2 = 1/r2;
        end

        if isnan(u3)
            u3r3 = 0;
            v3r3 = 0;
            nm3 = 0;
        else
            u3r3 = u3/r3;
            v3r3 = v3/r3;
            nm3 = 1/r3;
        end

        if isnan(u4)
            u4r4 = 0;
            v4r4 = 0;
            nm4 = 0;
        else
            u4r4 = u4/r4;
            v4r4 = v4/r4;
            nm4 = 1/r4;
        end
        norm = nm1 + nm2 + nm3 + nm4;
        norm = 1/norm;

        u0 = u1r1 + u2r2 + u3r3 + u4r4;
        v0 = v1r1 + v2r2 + v3r3 + v4r4;
        u0 = u0*norm;
        v0 = v0*norm;



%         norm = 1/r1+1/r2+1/r3+1/r4;
%         norm = 1/norm;
%         if isnan(u1) || isnan(u2) || isnan(u3) || isnan(u4)
%             break
%         end
%         u0 = u1/r1 + u2/r2 + u3/r3 + u4/r4;
%         v0 = v1/r1 + v2/r2 + v3/r3 + v4/r4;
%         u0 = u0*norm;
%         v0 = v0*norm;
    end

    U = sqrt(u0^2 + v0^2);
    ddx = abs(x2 - x1);
    dt = ddx/U;
    dt = dt*dr;

    startx = startx + dt*u0;
    starty = starty + dt*v0;

    X(i) = startx;
    Y(i) = starty;
%     plot(X(i),Y(i),'k.')
%     drawnow

%     toc
%     disp('b')
%     tic
    r1 = sqrt( (startx-x1).^2 + (starty-y1).^2 );
    r2 = sqrt( (startx-x2).^2 + (starty-y2).^2 );
    r3 = sqrt( (startx-x3).^2 + (starty-y3).^2 );
    r4 = sqrt( (startx-x4).^2 + (starty-y4).^2 );
    if r2 < r1 && r2 < r3 && r2 < r4
%     if j == 2
        row = row2;
    elseif r3 < r1 && r3 < r2 && r3 < r4
        col = col2;
    elseif r4 < r1 && r4 < r2 && r4 < r3
        row = row2;
        col = col2;
    end
      


end
