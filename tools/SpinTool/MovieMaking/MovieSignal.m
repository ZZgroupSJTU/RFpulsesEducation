function F = MovieSignal(signal,labelX,Nstep)

figHnd=figure
figPos = get(figHnd,'Position')
set(gcf,'color','white')

xlabel(labelX);

N=length(signal);

for k = 1:Nstep:N
    plot(signal(1:k),'k');
    axis([0 N min(signal) max(signal)]);
    set(gca,'Visible','off');
    F(k) = getframe(figHnd,[1 1 figPos(3) figPos(4)]);
end