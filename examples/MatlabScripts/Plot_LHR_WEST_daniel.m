close all

ne=WestHtoHD0.ne*20;
nHi=WestHtoHD0.nHi; nHi=0*nHi;
nDi=ne;%WestHtoHD0.nDi+nHi; 
nH2i=WestHtoHD0.nH2i;
nHDi=WestHtoHD0.nHDi;
nD2i=WestHtoHD0.nD2i+nHDi+nH2i;nD2i=nHi;nH2i=nHi;nHDi=nHi;

Rpos = WestHtoHD0.Rpos; %interp(Rpos,10);

R=250.0
a=67.0
b=70.0
lHFS=(R-184.3)/a
lLFS=(312.4-R)/a

rmin=(R-a*1.05); rmax=(R+a*1.05);

f = 48e6; % Hz
B0 = .38 * 1e4; % G
B=B0./Rpos*R;

HtoHD=0.01
radant=310

time=0.046

mycolor=lines(length(Rpos));

h=figure(1)
    set(h,'Units', 'Normalized', 'OuterPosition', [0.1-0.01 0.1-0.01 0.6 0.8]);
    Z=1;mu=1;
    wce=1.7588e7*B;%1.76
    wcHi=9.579e3*B*Z/1;%9.58
    wcDi=9.579e3*B*Z/2;
    wcH2i=9.579e3*B*Z/2;
    wcHDi=9.579e3*B*Z/3;
    wcD2i=9.579e3*B*Z/4;
    wcH3i=9.579e3*B*Z/3;
    wcHeII=9.579e3*B*Z/4;
    wcHeIII=9.579e3*B*2/4;

    Rpos_ = Rpos;%interp(Rpos,10);

    wpe=5.64e4*sqrt(ne);
    wpHi=1.32e3*Z*sqrt(nHi);
    wpDi=1.32e3*Z*sqrt(nDi/2);
    wpH2i=1.32e3*Z*sqrt(nH2i/2);
    wpHDi=1.32e3*Z*sqrt(nHDi/3);
    wpD2i=1.32e3*Z*sqrt(nD2i/4);

    epsperp=real(1 - ((wpe.^2)./((2*pi*f+i*1e5)^2-wce.^2)) ...
        - ((wpHi.^2)./((2*pi*f+i*1e4)^2-wcHi.^2)) ...
        - ((wpDi.^2)./((2*pi*f+i*1e4)^2-wcDi.^2)) ...
        - ((wpH2i.^2)./((2*pi*f+i*1e4)^2-wcH2i.^2)) ...
        - ((wpHDi.^2)./((2*pi*f+i*1e4)^2-wcHDi.^2)) ...
        - ((wpD2i.^2)./((2*pi*f+i*1e4)^2-wcD2i.^2)));
    locs = (find(diff(sign(epsperp))==2 | diff(sign(epsperp))==-2));

    ax1=subplot(3,1,1)
    y=ne;
    plot(Rpos,y,'.-','LineWidth',2,'Color',mycolor(1,:)); grid on; hold on
    xline(radant,'.-','Color',mycolor(3,:),'LineWidth',4);
    if (length(locs)>0)
        xline(Rpos_(locs(1)),'--','Color',mycolor(2,:),'LineWidth',2);
    end
    xline(R-lHFS*a,'k-');
    xline(R+lLFS*a,'k-');
    if (length(locs)>1)
        for j=2:length(locs)
            xline(Rpos_(locs(j)),'--','Color',mycolor(2,:),'LineWidth',2);
        end
    end
    title(['Tomator-1D/WEST 47MHz HtoHD=' num2str(HtoHD) ',  ne(t=' num2str(time*1e3,'%10.2f') 'ms)'])
    legend('ne(r)','Antenna','LHR','Limiter','Location','NorthWest')
    ylabel('ne [1/cm3]')
%         xlabel('R pos [cm]')
    % ylim([0 15e10])
    xlim([rmin rmax])
    set(gca,'FontSize',16) % right y-axis
    hold off

%     ax2=subplot(3,1,2)
%     y=Te(:,ind);
%     plot(Rpos,y,'.-','LineWidth',2,'Color',mycolor(1,:)); grid on; hold on
%     xline(radant,'.-','Color',mycolor(3,:),'LineWidth',4);
%     if (length(locs)>0)
%         xline(Rpos_(locs(1)),'--','Color',mycolor(2,:),'LineWidth',2);
%     end
%     xline(R-lHFS*a,'k-');
%     xline(R+lLFS*a,'k-');
%     if (length(locs)>1)
%         for j=2:length(locs)
%             xline(Rpos_(locs(j)),'--','Color',mycolor(2,:),'LineWidth',2);
%         end
%     end
%         title(['Tomator-1D/KIPT-RF TOMAS H2 25MHz Te(t=' num2str(time(ind)*1e3,'%10.2f') 'ms)'])
%     legend('Te(r)','Antenna','LHR','Limiter','Location','NorthWest')
%     ylabel('Te [eV]')
% %         xlabel('R pos [cm]')
%     ylim([0 6])
%     xlim([rmin rmax])
%     set(gca,'FontSize',16)
%     hold off

	ax3=subplot(3,1,3)
        plot(Rpos,epsperp,'LineWidth',2); grid on; hold on
        ylim([-1 1])

%     y=([PRFe(:,ind)  PRFHi(:,ind)  PRFH2i(:,ind)  PRFH3i(:,ind)  PRFHeII(:,ind)  PRFHeIII(:,ind)])*1.6e-19;
%     plot(Rpos,y,'LineWidth',2); grid on; hold on
%     xline(radant,'.-','Color',mycolor(3,:),'LineWidth',4);
%     if (length(locs)>0)
%         xline(Rpos_(locs(1)),'--','Color',mycolor(2,:),'LineWidth',2);
%     end
%         xline(Rpos_(locs(1)),'--','Color',mycolor(2,:),'LineWidth',2);
%     xline(R-lHFS*a,'k-');
%     xline(R+lLFS*a,'k-');
%     if (length(locs)>1)
%         for j=2:length(locs)
%             xline(Rpos_(locs(j)),'--','Color',mycolor(2,:),'LineWidth',2);
%         end
%     end
%         title('power deposition profiles(t_{end})')
%     legend('e','Hi','H2i','H3i','HeII','HeIII','Location','NorthWest')
%     ylabel('power [eV/s/cm3]')
%     ylim([0 0.04])
    xlim([rmin rmax])
%     xlabel('R pos [cm]')
%     set(gca,'FontSize',16)
%     hold off

    drawnow
%         % Capture the plot as an image
%         frame = getframe(h);
%         im = frame2im(frame);
%         if savevid
%             imwrite(im,[filename '_' num2str(ind)],'png');
%         end

    linkaxes([ax1 ax2 ax3],'x')
    
    