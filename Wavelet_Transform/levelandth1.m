function snrwave=levelandth1(y,s,wave,n)
for j=1:n
    xdh=wden(s,'sqtwolog','s','mln',5,wave(j,:));
    snrwave(j)=snrr(y,xdh);
end
end
