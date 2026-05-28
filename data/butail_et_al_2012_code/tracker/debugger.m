function debugger(Xh, Xi, k,  p, frmloc, imstream, fg, cams, tc, stage)
% function debugger(k, scan, p, frmloc, imstream, fg, cams, tc, Xi, Np, stage)
% nc=size(imstream,2);
% global gf;
% switch stage
%     case 'pre_update'
%         for cc=1:nc
%             figure(cc); gcf; clf;
%             imshow([frmloc, imstream(cc).flist(k).name]);
%             hold on;
%             figure(cc+2); gcf; clf;
% %             imshow(1-fg(:,:,cc)); 
%             imshow(1-fg{cc}); 
%             hold on;
%             for cl=1:size(scan(k-1).C,2)
%                 for t=scan(k-1).C(cl).t_id
%                     [tr br]=id2tbr(t,gf);
%                     [r c]=getind(Xi.nX, k, tr, 1:Xi.nX, Np);
%                     try
%                     tpzh=w2cam(p(r(Xi.ri),c, br), cams(cc));
%                     catch
%                         keyboard
%                     end
%                     tv=tra2b(p(r(Xi.rdi),c, br), [cams(cc).trm(1:3,1:3), [0 0 0]'; 0 0 0 1]);
%                     ntv=normv(tv);
%                     figure(cc); gcf;
%                     plot(tpzh(1,:), tpzh(2,:), '.', 'Color', tc(tr,:), 'MarkerSize', 2)
%                     quiver(tpzh(1,1:10:end), tpzh(2,1:10:end), ntv(1,1:10:end), ntv(2,1:10:end),...
%                                     'Color', tc(tr,:));
%                     figure(cc+2); gcf;
%                     
%                 end
%             end
%             figure(cc); gcf;
%             xlabel(sprintf('CAM %d PRE-Update ....', cc));
%             figure(cc+2); gcf;
%             xlabel(sprintf('CAM %d PRE-Update ....', cc));
%             drawnow;
%         end
%         keyboard;
%     case 'post_update'
%         for cc=1:nc
%             figure(cc); gcf;
%             imshow([frmloc, imstream(cc).flist(k).name]);
%             hold on;
%             for t_ii=ids'
%                 r=snipx(Xi.cX, t_ii);
%                 tpzh=w2cam(p(r(Xi.ri),:), cams(cc));
%                 tv=tra2b(p(r(Xi.rdi),:), [cams(cc).trm(1:3,1:3), [0 0 0]'; 0 0 0 1]);
%                 ntv=normv(tv);
%                 
%                 nzi=find(Xh(r(Xi.ri(1)),:)~=0);
%                 tzh=w2cam(Xh(r(Xi.ri),nzi(1):k), cams(cc));
%                 plot(tzh(1,:), tzh(2,:), 'Color', tc(t_ii,:), 'linewidth', 1);
%                 quiver(tpzh(1,1:10:end), tpzh(2,1:10:end), ntv(1,1:10:end), ntv(2,1:10:end), 'Color', tc(t_ii,:));
%             end
%             xlabel(sprintf('CAM %d POST-Update ....', cc), 'Color', 'r');
%             
%             figure(cc+2); gcf;
%             for t_ii=ids'
%                 r=snipx(Xi.cX, t_ii);
%                 nzi=find(Xh(r(Xi.ri(1)),:)~=0);
%                 tzh=w2cam(Xh(r(Xi.ri),nzi(1):k), cams(cc));
%                 tv=tra2b(Xh(r(Xi.rdi),k), [cams(cc).trm(1:3,1:3), [0 0 0]'; 0 0 0 1]);
%                 ntv=normv(tv);
%                 plot(tzh(1,:), tzh(2,:), '-o', 'Color', tc(t_ii,:),...
%                             'MarkerSize', 10, 'LineWidth', 2);
%                 quiver(tzh(1,end), tzh(2,end), ntv(1), ntv(2), 30, 'Color', tc(t_ii,:));
%             end
%             
%             % show the ones that were terminated in this time-step
%             idprev=getids(Xh,k-1,Xi);
%             idterminate=idprev(~ismember(idprev,ids));
%             for t_ii=idterminate'
%                 r=snipx(Xi.cX, t_ii);
%                 nzi=find(Xh(r(Xi.ri(1)),:)~=0);
%                 tzh=w2cam(Xh(r(Xi.ri),nzi(1):k-1), cams(cc));
%                 tv=tra2b(Xh(r(Xi.rdi),k-1), [cams(cc).trm(1:3,1:3), [0 0 0]'; 0 0 0 1]);
%                 ntv=normv(tv);
%                 plot(tzh(1,:), tzh(2,:), '-o', 'color', ones(1,3)*.5, ...
%                             'MarkerSize', 5, 'LineWidth', 1);
%                 quiver(tzh(1,end), tzh(2,end), ntv(1), ntv(2), 30, 'color', ones(1,3)*.5);
%             end
%             
%             xlabel(sprintf('CAM %d POST-Update ....', cc), 'Color', 'r');
%             drawnow;
%         end
%         keyboard;
%     case 'new_targets'
%         ids=getids(Xh,k+1,Xi);
%         for cc=1:nc
%             figure(cc); gcf; clf;
%             imshow([frmloc, imstream(cc).flist(k).name]);
%             hold on;
%             for t_ii=ids'
%                 r=snipx(Xi.cX, t_ii);
%                 tpzh=w2cam(p(r(Xi.ri),:), cams(cc));
%                 nzi=find(Xh(r(Xi.ri(1)),:)~=0);
%                 tzh=w2cam(Xh(r(Xi.ri),nzi(1):k+1), cams(cc));
%                 tv=tra2b(p(r(Xi.rdi),:), [cams(cc).trm(1:3,1:3), [0 0 0]'; 0 0 0 1]);
%                 ntv=normv(tv);
%                 plot(tzh(1,:), tzh(2,:), 'Color', tc(t_ii,:), 'linewidth', 1);
%                 plot(tpzh(1,:), tpzh(2,:), '.', 'Color', tc(t_ii,:), 'MarkerSize', 2)
%                 quiver(tpzh(1,1:10:end), tpzh(2,1:10:end), ntv(1,1:10:end), ntv(2,1:10:end), 'Color', tc(t_ii,:));
%             end
%             xlabel(sprintf('CAM %d New-targets ....', cc));
% 
%             figure(cc+2); gcf; clf;
% %             imshow(1-fg(:,:,cc)); hold on;
%             imshow(1-fg{cc}); hold on;
%             for t_ii=ids'
%                 r=snipx(Xi.cX, t_ii);
%                 tzh=w2cam(Xh(r(Xi.ri),k+1), cams(cc));
%                 tv=tra2b(Xh(r(Xi.rdi),k+1), [cams(cc).trm(1:3,1:3), [0 0 0]'; 0 0 0 1]);
%                 ntv=normv(tv);
%                 plot(tzh(1), tzh(2), '+', 'Color', tc(t_ii,:)+.75*(ones(1,3)-tc(t_ii,:)),...
%                             'MarkerSize', 10, 'LineWidth', 2);
%                 quiver(tzh(1), tzh(2), ntv(1), ntv(2), 30, 'Color', tc(t_ii,:)+.75*(ones(1,3)-tc(t_ii,:)));
% 
%                 text(tzh(1)+1, tzh(2)+1, sprintf('%d', t_ii), 'Color', tc(t_ii,:)+.75*(ones(1,3)-tc(t_ii,:)));
%             end
%             xlabel(sprintf('CAM %d NEW-targets ....', cc));
%             drawnow;
%         end
%         keyboard;        
% end