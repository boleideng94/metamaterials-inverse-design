function ang_b = ang_bet(u,v)
%     u = [u;0];
%     v = [v;0];
%     ang_b = atan2(norm(cross(u,v)),dot(u,v));
    CosTheta = (dot(u,v))/(norm(u)*norm(v));
    ang_b = acos(CosTheta);
end