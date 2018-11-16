function [ind] = plotcsensorth( Sens_orth,params, removed , lb ,Np)
    %Take matrix A where A_ij = j if parameter i is above the threshold for
    %experiment j
    rem = remaining(Np , removed);
    Np=Np - length(removed);
    fs=12;
    figure(); clf;
    %set(gcf,'Position',[234     3   564   695-105],'PaperPositionMode','auto');
    %[pl,ind]=sort(sum((A~=0)'))
    [pl,ind]=sort(mean(log10(Sens_orth'))); % ind dicates the ranking of the parameters based on sens. magnitude
    barh([(pl)])
    set(gca,'YTick',1:Np);
    set(gca,'Position',[0.2 0.1 0.75 0.85])
    %set(gca, 'YTickLabel', params);
    [hx,hy] = format_ticks(gca,' ',params(rem(ind)),lb:1:0,[],0,0,0.05,'FontSize',fs,'FontWeight','Bold');
    %format_ticks(gca,' ',params(ind),0:400:1600,[],0,0,0.05,'FontSize',fs,'FontWeight','Bold');
    set(gca,'FontSize',fs);
    xlabel('sensitivity magnitude','FontSize',fs)
    %xlabel('Cardinality','FontSize',fs)
end

