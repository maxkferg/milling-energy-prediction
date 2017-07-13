function plot_variance(x,lower,upper,color)
    set(fill([x,x(end:-1:1)],[upper,fliplr(lower)],color),'EdgeColor','none');
    plot(x,upper,color);
    plot(x,lower,color);
end