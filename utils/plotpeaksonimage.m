function plotpeaksonimage(peakim,im,pflag,color,ysrf,ybtm)

figure; axis ij; axis on; hold on;
[y,x] = find(peakim > 0);
if pflag == 'p'
    imagesc(im);colormap('gray');
    plot(x,y,'.','color',color,'markersize',1);
else
ind = sub2ind(size(peakim),y,x);
val = peakim(ind);%max(peakim(:)) - peakim(ind) + 1;
scatter(x,y,10,val,'.');
caxis([0 150]);colormap('cool');
plot(1:length(ysrf),ysrf,'-k','linewidth',3);
plot(1:length(ybtm),ybtm,'-k','linewidth',3);
colorbar;%caxis([0 max(peakim(:))]);
end
ylim([0 1200]); xlim([0 size(im,2)]);

plot([100 100],[30 100],'k-','linewidth',3); % 200 m
plot([100 478],[30 30],'k-','linewidth',3); % 500 m
set(findall(gcf,'-property','FontSize'),'FontSize',10);
end