function s = Area3D4P(A,B,C,D)
s1 = helen(A,B,C);
s2 = helen(A,C,D);
s3 = helen(A,B,D);
s4 = helen(B,C,D);
s = (s1+s2+s3+s4)/2;