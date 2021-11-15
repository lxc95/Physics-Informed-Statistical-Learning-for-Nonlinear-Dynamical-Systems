function s = helen(x,y,z)
a = lenth(x,y);
b = lenth(x,z);
c = lenth(y,z);
p = (a+b+c)/2;
s = sqrt(p*(p-a)*(p-b)*(p-c));