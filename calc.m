function [prec,rec,acc] = calc(CM)
a = CM(1,1);
b = CM(1,2);
c = CM(2,1);
d = CM(2,2);

if(a+d > b+c)
   prec = ((a/(a+c)) + (d/(b+d)))/2;
   rec = ((a/(a+b)) + (d/(d+c)))/2;
   acc = (a+d)/(a+b+c+d);
else
   prec = ((b/(b+d)) + (c/(a+c)))/2;
   rec = ((b/(a+b)) + (c/(c+d)))/2;
   acc = (b+c)/(a+b+c+d);
end