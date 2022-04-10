function [z,snrxd]=level(y,s,th,wave)%sÎª¼ÓÔëĞÅºÅ
for k=1:10
xdh=wden(s,th,'s','mln',k,wave);
snrxd(k)=snrr(y,xdh);
end
[~,z]=max(snrxd);
end
