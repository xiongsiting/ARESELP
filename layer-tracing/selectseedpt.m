function seedpt = selectseedpt(peakim)
a = peakim(peakim > 0);
[x,y] = find(peakim > 0);
indices = sub2ind(size(peakim),x,y);
val = peakim(indices);
seedpt = [val,x,y];
seedpt = flipud(sortrows(seedpt));
pd = fitdist(a,'lognormal');
thresh = exp(pd.mu + pd.sigma.^2 /2);
seedpt(seedpt(:,1) < thresh,:) = [];
seedpt(:,1) = [];
end
