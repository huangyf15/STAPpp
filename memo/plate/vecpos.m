function res=vecpos(a,b)
if (b>=a)
    res=b*(b+1)/2-a;
else
    res=a*(a+1)/2-b;
end