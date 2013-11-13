function costM = calculateDTW(signals)
    
    tic;
    fid = fopen('distances/costs.data','w');
    fprintf(fid,datestr(now));
    fprintf(fid,'\n\nMFCC Spectrogram\n');
    
    gcf=figure('Color',[1 1 1]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [0.25 2.5 8 6]);
       
    costM = zeros(length(signals),length(signals));
    
    for i=1:length(signals)
        for j=i:length(signals)
            signal1 = signals(i).wav;
            signal2 = signals(j).wav;
%             display('Computing similarity matrix');
            M = simmx(abs(signal1),abs(signal2));
%             display('Dinamic programming (DTW)');
            [p q C] = dp2(1-M);
            
            subplot(121);
            imagesc(M);
            colormap(1-gray);
            pbaspect([1 1 1]);

            hold on;
            plot(p,q,'r');
            hold off;

            subplot(122);
            imagesc(C);
            pbaspect([1 1 1]);
            hold on;
            plot(p,q,'r');
            hold off;
            
            costM(i,j) = C(size(C,1), size(C,2));
            costM(j,i) = C(size(C,1), size(C,2));
            
            fprintf(fid,['Cost from ' num2str(i) ' to ' num2str(j) ...
                                ': ' num2str(C(size(C,1), size(C,2))) '\n']);
            saveas(gcf,['distances/MFCC_' num2str(i) '_to_' num2str(j) '.png'],'png');
        end
    end
    
    display(costM);
    toc
end