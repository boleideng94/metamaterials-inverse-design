function centroid_p = hole_cad(fid,p1,p2,p3,p4,d,flag)

polyin = polyshape([p1(1),p2(1),p3(1),p4(1)],[p1(2),p2(2),p3(2),p4(2)]);
[x,y] = centroid(polyin);
centroid_p = [x,y];

% centroid_p = 

vec1 = (p2-p1);
vec2 = (p4-p1);
p11 = p1 + vec1.*(-d/2/vec1(1));
p12 = p1 + vec2.*(-d/2/vec2(1));

vec1 = (p1-p2);
vec2 = (p3-p2);
p21 = p2 + vec1.*(-d/2/vec1(2));
p22 = p2 + vec2.*(-d/2/vec2(2));

cross0 = cross([-vec1,0],[vec2,0]);

if cross0(3)<0 || abs(vec1(2)*vec2(2))<0.05
    p21 = p2 + d/2.*[1,-1];
    p22 = p2 + d/2.*[-1,-1];    
end

vec1 = (p2-p3);
vec2 = (p4-p3);
p31 = p3 + vec1.*(d/2/vec1(1));
p32 = p3 + vec2.*(d/2/vec2(1));

vec1 = (p1-p4);
vec2 = (p3-p4);
p41 = p4 + vec1.*(d/2/vec1(2));
p42 = p4 + vec2.*(d/2/vec2(2));

if p31(2)>p22(2)
    p31(2) = p21(2);
    p22(1) = p31(1);
end

if p4(2)>p3(2)
    p42 = p4 + [-d/2,d/2];
end

if p4(2)>p1(2)
     p41 = p4 + [d/2,d/2];
end   

% if p42(2)>p32(2)
    

if flag == 0
    fprintf(fid,linecad(p31(1),p31(2),p32(1),p32(2)));
    fprintf(fid,linecad(p22(1),p22(2),p31(1),p31(2)));
    fprintf(fid,linecad(p32(1),p32(2),p42(1),p42(2)));
    fprintf(fid,linecad(p11(1),p11(2),p12(1),p12(2)));
    fprintf(fid,linecad(p11(1),p11(2),p21(1),p21(2)));
    fprintf(fid,linecad(p12(1),p12(2),p41(1),p41(2)));
    fprintf(fid,linecad(p21(1),p21(2),p22(1),p22(2)));
    fprintf(fid,linecad(p41(1),p41(2),p42(1),p42(2)));    
end

if flag == 1
    fprintf(fid,linecad(p11(1),p11(2),p12(1),p12(2)));
    fprintf(fid,linecad(p11(1),p11(2),p2(1),p2(2)));
    fprintf(fid,linecad(p12(1),p12(2),p4(1),p4(2)));
end

if flag == 2
    fprintf(fid,linecad(p31(1),p31(2),p32(1),p32(2)));
    fprintf(fid,linecad(p2(1),p2(2),p31(1),p31(2)));
    fprintf(fid,linecad(p32(1),p32(2),p4(1),p4(2)));   
end
end