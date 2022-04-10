function  [Xs,Ps] = RTS(F,xf,Pf,xf1,Pf1)
    Xs = xf1;Ps = Pf1;
    for k = length(F)-1:-1:1
        Ks = Pf1{k}*F{k}'*invbc(Pf{k+1});
        Xs{k} = xf1{k} + Ks*(Xs{k+1}-xf{k+1});
        Ps{k} = Pf1{k} + Ks*(Ps{k+1} - Pf{k+1})*Ks';
    end
end