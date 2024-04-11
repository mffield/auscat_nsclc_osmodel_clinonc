function plotCumulHist(model,res_str, pi_str, col)


    samp = 200;
    xi = linspace(0,1,samp);
    yi = []; tyi=[];

    for i = 1:length(model.client); client_y(i).yi=[]; end
    for k = 1:length(model.client(1).(res_str))
        histlvl=0; cumhist=0;

        for i = 1:length(model.client)
            if length(model.client(i).(res_str)(k).(pi_str).CumHist)>1
                CH = [model.client(i).(res_str)(k).(pi_str).CumHist;]; 
                CHL = [model.client(i).(res_str)(k).(pi_str).CumHistLevels; ]; %model.client(i).validation_result(k).(pi_str).minmax(2)+eps
                CHN = model.client(i).(res_str)(k).(pi_str).NumSamples;
                cumhist_c = [0; CHN*model.client(i).(res_str)(k).(pi_str).minmax(1);
                    diff(CH*CHN)];

                for j=1:length(CHL)
        %             if cumhist_c(j)>0
                        histlvl = sort([histlvl CHL(j)]); ind=find(histlvl==CHL(j),1,'first');
                        cumhist = [cumhist(1:ind-1) cumhist_c(j) cumhist(ind:end)];
        %             end
                end
                if length(CHL)>1
                    client_y(i).yi = [client_y(i).yi; interp1(CHL,double((cumsum(CH)/sum(CH))),xi,'linear','extrap')];    
                    client_y(i).yi(client_y(i).yi>1)=1;
                    client_y(i).yi(client_y(i).yi<0)=0;
                end
            end
        end
        if length(histlvl)>1
            tyi = [tyi; interp1(histlvl,double((cumsum(cumhist)/sum(cumhist))),xi,'linear','extrap')];
            tyi(tyi>1)=1;
        end
    end

    ptyi = prctile(tyi,[2.5,50,97.5]);
    plot(xi,ptyi(2,:),'Linewidth',2,'Color',col);
    hold on; 
    patch([xi'; flipud(xi')], [ptyi(3,:)'; flipud(ptyi(1,:)')], col, 'FaceAlpha',0.2, 'EdgeColor','none','HandleVisibility','off')


end