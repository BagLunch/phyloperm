ammo = unique(Mcgowan2004data,'rows');
tr1 = ammo.D;
tr2 = ammo.W;
[tx,ty] = minboundtri(tr1,tr2);
hold off
plot(tr1,tr2,'ro',tx,ty,'b-')
hold on
[K,A] = convhull(tr1,tr2);
plot(tr1(K),tr2(K))
A/polyarea(tx,ty)

ratios = [];
for i = 1:1000
    tr1 = tr1(randperm(length(tr1)));
    tr2 = tr2(randperm(length(tr2)));
    [tx,ty] = minboundtri(tr1,tr2);
    [K,A] = convhull(tr1,tr2);
    ratios = [ratios,A/polyarea(tx,ty)];
end
csvwrite("ratios_operm.csv",ratios)

sum(tr2 > 1./tr1)
counts = [];
for i = 1:1000
    tr1 = tr1(randperm(length(tr1)));
    tr2 = tr2(randperm(length(tr2)));
    counts = [counts,sum(tr2 > 1./tr1)];
end

%phyloperms
ammoDphyloperms = table2array(ammoDphyloperms);
ammoWphyloperms = table2array(ammoWphyloperms);
ratios = [];
for i = 1:1000
    tr1 = ammoDphyloperms(:,i);
    tr2 = ammoWphyloperms(:,i);
    [tx,ty] = minboundtri(tr1,tr2);
    [K,A] = convhull(tr1,tr2);
    ratios = [ratios,A/polyarea(tx,ty)];
end
csvwrite("ratios_phyloperm.csv",ratios)

