function re=linecad(cx1,cy1,cx2,cy2)
re=sprintf('(command "line" "%g,%g" "%g,%g") \n',cx1,cy1,cx2,cy2);
% re=strcat('(command "circle" ','"',num2str(cx1),',',num2str(cy1),'"\n','"',num2str(cx2),',',num2str(cy2),'"',' ',')');
end