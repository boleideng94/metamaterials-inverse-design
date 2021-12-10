clear all
a = 10;
d = 1;
bt = 5;
mt = 5;
mout = 5;
TOL = 0.06;
hole_r = 1.3;
hole_r2 = 0.5;

rc_m = 2; rc_n = 2;
p_m = 5; p_n = 4;
m = rc_m*p_m;
n = rc_m*p_n+1;

sample_height = n*a+2*bt;
sample_width = m*a;

case_name = 'Fig1_random';
input = load(['input_design/',case_name,'.txt']);

rc_rand = reshape(input,[4,rc_m,rc_n]);

rc_rand(:,1:rc_n,rc_m+1) = rc_rand(:,1:rc_n,1);
rc_rand(:,rc_n+1,1:rc_m) = rc_rand(:,1,1:rc_m);
n = rc_n*p_n+2; m = rc_m*p_m;

hori_shift =zeros([2,n,m+1]);
vert_shift = zeros([2,n+1,m]);

for p_i = 1:p_n
    for p_j = 1:p_m
        vert_shift(:,(p_i-1)*rc_n+2:(p_i)*rc_n+2,(p_j-1)*rc_n+1:(p_j)*rc_n+1) = rc_rand(1:2,:,:);
        hori_shift(:,(p_i-1)*rc_n+2:(p_i)*rc_n+2,(p_j-1)*rc_n+1:(p_j)*rc_n+1) = rc_rand(3:4,:,:);
    end
end
hori_shift(:,n,:) = 0;

name = ['autocad_design/',case_name,'.scr'];
fid = fopen(name,'w');

comd = strcat("(command ","setvar ",'"osmode"', " 0)\n");
fprintf(fid,sprintf(comd)); 

% for pillar holes
for i=1:n-1
    for j=1:m-1
        center = [j,i].*a;
            p1 = center + ([0.5,0] + vert_shift(:,i+1,j+1)').*(a+2*TOL);
            p2 = center + ([0,0.5] + hori_shift(:,i+1,j+1)').*(a+2*TOL);
            p3 = center + ([-0.5,0] + vert_shift(:,i+1,j)').*(a+2*TOL);
            p4 = center + ([0,-0.5] + hori_shift(:,i,j+1)').*(a+2*TOL);
            centroid_p = hole_cad(fid,p1,p2,p3,p4,d,0);
            fprintf(fid,circle_cad(centroid_p(1),centroid_p(2),hole_r));
    end
end


% for other parts
j = 0;
for i=1:n-1
        center = [j,i].*a;
            p1 = center + ([0.5,0] + vert_shift(:,i+1,j+1)').*(a+2*TOL);
            p2 = center + ([0,0.5] + hori_shift(:,i+1,j+1)').*(a) + [TOL,0];
            p3 = center + ([-0.5,0] + [0,0]).*(a+2*TOL);
            p4 = center + ([0,-0.5] + hori_shift(:,i,j+1)').*(a) + [TOL,0];
            hole_cad(fid,p1,p2,p3,p4,d,1);
end

j = m;
for i=1:n-1
        center = [j,i].*a;
            p1 = center + ([0.5,0] + [0,0]).*(a+2*TOL);
            p2 = center + ([0,0.5] + hori_shift(:,i+1,j+1)').*(a)+ [-TOL,0];
            p3 = center + ([-0.5,0] + vert_shift(:,i+1,j)').*(a+2*TOL);
            p4 = center + ([0,-0.5] + hori_shift(:,i,j+1)').*(a)+ [-TOL,0];
            hole_cad(fid,p1,p2,p3,p4,d,2);
end

b1 = [0,n-0.5].*a + [TOL,0];
b2 = b1 + [0,bt-TOL];
b3 = [m,n-0.5].*a + [-TOL,0];
b4 = b3 + [0,bt-TOL];
fprintf(fid,linecad(b1(1),b1(2),b2(1),b2(2)));
fprintf(fid,linecad(b4(1),b4(2),b2(1),b2(2)));
fprintf(fid,linecad(b3(1),b3(2),b4(1),b4(2)));
b1 = [0,0.5].*a + [TOL,0];
b2 = b1 - [0,bt-TOL];
b3 = [m,0.5].*a + [-TOL,0];
b4 = b3 - [0,bt-TOL];
fprintf(fid,linecad(b1(1),b1(2),b2(1),b2(2)));
fprintf(fid,linecad(b4(1),b4(2),b2(1),b2(2)));
fprintf(fid,linecad(b3(1),b3(2),b4(1),b4(2)));

inner_box = [-mt-TOL,a/2-bt-mt-TOL;sample_width + mt+TOL,sample_height+a/2-bt+mt+TOL];
rect_cad(fid,inner_box(1,:),inner_box(2,:));

center = [inner_box(1,1) + hole_r*2,(inner_box(1,2)+inner_box(2,2))/2];
fprintf(fid,circle_cad(center(1),center(2),hole_r));
center = [inner_box(2,1) - hole_r*2,(inner_box(1,2)+inner_box(2,2))/2];
fprintf(fid,circle_cad(center(1),center(2),hole_r));
center = [(inner_box(1,1)+inner_box(2,1))/2,inner_box(1,2) + hole_r*2];
fprintf(fid,circle_cad(center(1),center(2),hole_r));
center = [(inner_box(1,1)+inner_box(2,1))/2,inner_box(2,2) - hole_r*2];
fprintf(fid,circle_cad(center(1),center(2),hole_r));

rect_cad(fid,[-mt-mout,a/2-bt-mt-mout],[sample_width + mt+mout,sample_height+a/2-bt+mt+mout]);

%  walls
shift = [sample_width+30,0];
j = 0;
for i=1:n-1
        center = [j,i].*a+shift;
            p1 = center + ([0.5,0] + vert_shift(:,i+1,j+1)').*a;
            p2 = center + ([0,0.5] + hori_shift(:,i+1,j+1)').*a;
            p3 = center + ([-0.5,0] + [0,0]).*a;
            p4 = center + ([0,-0.5] + hori_shift(:,i,j+1)').*a;
            hole_cad(fid,p1,p2,p3,p4,d,1);
end
j = m;
for i=1:n-1
        center = [j,i].*a+shift;
            p1 = center + ([0.5,0] + [0,0]).*a;
            p2 = center + ([0,0.5] + hori_shift(:,i+1,j+1)').*a;
            p3 = center + ([-0.5,0] + vert_shift(:,i+1,j)').*a;
            p4 = center + ([0,-0.5] + hori_shift(:,i,j+1)').*a;
            hole_cad(fid,p1,p2,p3,p4,d,2);
end


b1 = [0,0.5].*a+shift;
b2 = b1 - [0,bt];
b3 = [0,n-0.5].*a+shift;
b4 = b3 + [0,bt];
b5 = b2 + [-mt,0];
b6 = b4 + [-mt,0]; 

fprintf(fid,linecad(b1(1),b1(2),b2(1),b2(2)));
fprintf(fid,linecad(b3(1),b3(2),b4(1),b4(2)));
fprintf(fid,linecad(b5(1),b5(2),b6(1),b6(2)));
fprintf(fid,linecad(b4(1),b4(2),b6(1),b6(2)));
fprintf(fid,linecad(b2(1),b2(2),b5(1),b5(2)));

st = b6+[0,10];
ed = st + [sample_width+2*mt,mt];
rect_cad(fid,st,ed);
st = b5-[0,10+mt];
ed = st + [sample_width+2*mt,mt];
rect_cad(fid,st,ed);

b1 = [n,0.5].*a+shift;
b2 = b1 - [0,bt];
b3 = [n,n-0.5].*a+shift;
b4 = b3 + [0,bt];
b5 = b2 + [mt,0];
b6 = b4 + [mt,0]; 
fprintf(fid,linecad(b1(1),b1(2),b2(1),b2(2)));
fprintf(fid,linecad(b3(1),b3(2),b4(1),b4(2)));
fprintf(fid,linecad(b5(1),b5(2),b6(1),b6(2)));
fprintf(fid,linecad(b4(1),b4(2),b6(1),b6(2)));
fprintf(fid,linecad(b2(1),b2(2),b5(1),b5(2)));

% rect_cad(fid,botleft,topright);

% pillars
shift = [2*sample_width+40,0];
for i=1:n-1
    for j=1:m-1
        center = [j,i].*a + shift;
            p1 = center + ([0.5,0] + vert_shift(:,i+1,j+1)').*a;
            p2 = center + ([0,0.5] + hori_shift(:,i+1,j+1)').*a;
            p3 = center + ([-0.5,0] + vert_shift(:,i+1,j)').*a;
            p4 = center + ([0,-0.5] + hori_shift(:,i,j+1)').*a;
            centroid_p = hole_cad(fid,p1,p2,p3,p4,d,0);
            fprintf(fid,circle_cad(centroid_p(1),centroid_p(2),hole_r2));
    end
end
