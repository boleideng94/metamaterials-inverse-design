function rect_cad(fid,botleft,topright)
fprintf(fid,linecad(botleft(1),botleft(2),topright(1),botleft(2))); 
fprintf(fid,linecad(botleft(1),botleft(2),botleft(1),topright(2))); 
fprintf(fid,linecad(botleft(1),topright(2),topright(1),topright(2))); 
fprintf(fid,linecad(topright(1),botleft(2),topright(1),topright(2))); 
end