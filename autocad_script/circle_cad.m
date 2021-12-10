function o = circle_cad(x1,y1,r)
    o = sprintf('(command "circle"  "%g,%g" "%g")\n', x1,y1,r);
end