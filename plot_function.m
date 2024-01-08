function [ ] = plot_function( x, y, xname,xsymbol,xunit, yname,ysymbol,yunit,legend_name,mode)

legend_entries = string(zeros(size(y,2),1));
h = find(mode);

if size(y,3) == 1
    plot(x,y)
elseif size(y,3) > 1
    plot(x(:,1),y(:,:,1));
    hold on;
    plot(x(:,1),y(:,:,2),'--');
end

title([yname,' ',ysymbol])
ylabel([yname,' ',ysymbol,' ',yunit])
xlabel([xname,' ',xsymbol,' ',xunit])

if size(y,2) > 1
    if size(y,3) == 1
        for k = 1:size(y,2)
%            legend_entries(k,1) = sprintf('%s%dM%d',legend_name,k,h(1));
        end
    elseif size(y,3) > 1
        for k = 1:size(y,2)*2
             if k <= size(y,2)
%                legend_entries(k,1) = sprintf('%s%dM%d',legend_name,k,h(1));
             elseif k > size(y,2)
%                 legend_entries(k,1) = sprintf('%s%dM%d',legend_name,k-size(y,2),h(2));
             end
        end
    end
else
%    legend_entries(1,1) = sprintf('%sM%d',legend_name,h(1));
end

%legend(legend_entries,'Location','NorthWest')



    
