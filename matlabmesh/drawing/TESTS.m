p1 = [2,3];
p2 = [4,5];
r = 3.4;
centers = circle2_2pr(p1, p2, r);

lineO = [2.5,3.5];
lineD = normalize([-1.2,1]);
hits = isect_line2_circle2(lineO, lineO+lineD, centers(1,:), r);

newplot;
hold all;
drawpoint(p1, 0.1);
drawpoint(p2, 0.1);
drawcircle( centers(1,:), r );
drawcircle( centers(2,:), r );

drawpoint(hits(1,:), 0.2);
drawpoint(hits(2,:), 0.2);
drawline(lineO-10*lineD,lineO+10*lineD);

axis equal;
hold off;


