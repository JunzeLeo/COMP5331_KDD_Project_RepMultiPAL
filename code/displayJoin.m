function [] = displayJoin(multi_joins, data, is_normalized)
    sublens = multi_joins.Ls;
    Cs = multi_joins.Cs;
    ids = multi_joins.ids;
    subseqs = multi_joins.subseqs;
    n_ids = size(ids,1);
    K = length(sublens);
    shift = 0;
    colors = ["#ADFF2F";"#FFE812";"#FF00FF";"#87CEFA";"#00FFFF";"#FF8383";"#FFFF00"];
    darkcolors = ["#006400";"#CD853F";"#8B008B";"#0000CD";"#008B8B";"#FF0000";"#DAA520"];
    M = size(data,2);
    XX = 1:M;
    legend_labels = strings(K,1);
    lgd = zeros(1, K);
    n_fig = 50;
    for i = 1:n_ids
        if mod(i,n_fig) == 1
            if i ~= 1
                %break;
                hold off;
                xlabel('Timeline');
                if is_normalized == 1
                    ylabel('Z-Normalized value (Shifted)');
                else
                    ylabel('Actual value (Shifted)');
                end
                figname = strcat(int2str(K),' Multi-way Joins, # of Joined Timeseries=',int2str(n_ids));
                title(figname);
                legend(lgd, cellstr(legend_labels));
                pause;
                close all;
            end
            shift = 0;
        end
        
        id = ids(i);
        if is_normalized == 1
            data(id,:) = zNorm(data(id,:));
        end
        if i ~= 1
            shift = shift+abs(min(0,min(data(id,:))));
        end
        mx = abs(max(data(id,:)));
        data(id,:) = shift+data(id,:);
        plot(XX,data(id,:),'Color','#696969');
        hold on;
        for j = 1:K
            subid = subseqs(i,j);
            if mod(i,n_fig) == 1 || j > 0
                legend_label = strcat('sublen=',int2str(sublens(j)),',C=',...
                num2str(Cs(j)));
                legend_labels(j) = legend_label;
                lgd(j) = plot(XX(subid:subid+sublens(j)-1),data(id,subid:subid+sublens(j)-1),...
                    'Color',colors(j),'LineWidth',8);
                plot(XX(subid:subid+sublens(j)-1),data(id,subid:subid+sublens(j)-1),...
                    'Color',darkcolors(j),'LineWidth',2);
            else
                plot(XX(subid:subid+sublens(j)-1),data(id,subid:subid+sublens(j)-1),...
                    'Color',colors(j),'LineWidth',8);
                 plot(XX(subid:subid+sublens(j)-1),data(id,subid:subid+sublens(j)-1),...
                    'Color',darkcolors(j),'LineWidth',2);
            end
            hold on;
        end
        shift = shift+mx+1.0;
    end
    hold off;
    xlabel('Timeline');
    if is_normalized == 1
        ylabel('Z-Normalized value (Shifted)');
    else
        ylabel('Actual value (Shifted)');
    end
    figname = strcat(int2str(K),' Multi-way Joins, # of Joined Timeseries=',int2str(n_ids));
    title(figname);
    legend(lgd, cellstr(legend_labels));
end
