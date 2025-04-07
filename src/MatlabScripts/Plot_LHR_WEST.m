clear all
close all
clc

file = 'Data/WEST/Res_20231127_103434.csv'

%%

opts = delimitedTextImportOptions("NumVariables", 61);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["tmain", "RadialPositions", "ne", "Ee", "dne", "dEe", "nue", "nH", "EH", "dnH", "dEH", "nuH", "nH2", "EH2", "dnH2", "dEH2", "nuH2", "nHi", "EHi", "dnHi", "dEHi", "nuHi", "nH2i", "EH2i", "dnH2i", "dEH2i", "nuH2i", "nH3i", "EH3i", "dnH3i", "dEH3i", "nuH3i", "nHeI", "EHeI", "dnHeI", "dEHeI", "nuHeI", "nHeII", "EHeII", "dnHeII", "dEHeII", "nuHeII", "nHeIII", "EHeIII", "dnHeIII", "dEHeIII", "nuHeIII", "nCI", "xnCI", "nCII", "xnCII", "nCIII", "xnCIII", "nCIV", "xnCIV", "nCV", "xnCV", "PRFe", "PRFHi", "PRFH2i", "PRFH3i", "PRFHeII", "PRFHeIII", "tnew"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
Res = readtable(['/home/ITER/wautert/Desktop/Tomator/tomator/' file], opts);
Data=Res;

%% Clear temporary variables
clear opts

R=250.0
a=67.0
b=70.0
lHFS=(R-184.3)/a
lLFS=(312.4-R)/a
radant=310
HtoHD=0.02

video=1; savevid=1;
every=10; % plot for every ... outputs
tmaxplot=1; % ignore
rmin=(R-a*1.05); rmax=(R+a*1.05);

kb=1.38e-23;
f = 48e6; % Hz
B0 = 3.8 * 1e4; % G
% ny=1;
% nz=-5:2:5;

pdp=1;  disper=1;

%% read data

m=301;n=length(Data.tmain)/m;

time    = reshape(Data.tmain,m,n);time=time(1,:);
Rpos    = reshape(Data.RadialPositions,m,n); Rpos=Rpos(:,1);
ne      = reshape(Data.ne,m,n);
Ee      = reshape(Data.Ee,m,n); Te=Ee./(3/2*ne);
% dne     = reshape(Data.(:,4),m,n);
% dEe     = reshape(Data.(:,5),m,n);
% nue    	= reshape(Data.(:,6),m,n);
nH     	= reshape(Data.nH,m,n);
EH     	= reshape(Data.EH,m,n); TH=EH./(3/2*nH);
% dnH     = reshape(Data.(,9),m,n);
% dEH     = reshape(Data.(,10),m,n);
% nuH    	= reshape(Data.(,11),m,n);
nH2    	= reshape(Data.nH2,m,n);
EH2     = reshape(Data.EH2,m,n); TH2=EH2./(3/2*nH2);
% dnH2    = reshape(Data.(,14),m,n);
% dEH2    = reshape(Data.(,15),m,n);
% nuH2    = reshape(Data.(,16),m,n);
nHi     = reshape(Data.nHi,m,n);
EHi     = reshape(Data.EHi,m,n); THi=EHi./(3/2*nHi);
% dnHi    = reshape(Data.(,19),m,n);
% dEHi    = reshape(Data.(,20),m,n);
% nuHi    = reshape(Data.(,21),m,n);
nH2i    = reshape(Data.nH2i,m,n);
EH2i    = reshape(Data.EH2i,m,n); TH2i=EH2i./(3/2*nH2i);
% dnH2i   = reshape(Data.(,24),m,n);
% dEH2i   = reshape(Data.(,25),m,n);
% nuH2i   = reshape(Data.(,26),m,n);
nH3i    = reshape(Data.nH3i,m,n);
EH3i    = reshape(Data.EH3i,m,n); TH3i=EH3i./(3/2*nH3i);
% dnH3i   = reshape(Data.(,29),m,n);
% dEH3i   = reshape(Data.(,30),m,n);
% nuH3i   = reshape(Data.(,31),m,n);
nHeI    = reshape(Data.nHeI,m,n);
EHeI    = reshape(Data.EHeI,m,n); THeI=EHeI./(3/2*nHeI);
% dnHeI   = reshape(Data.(,34),m,n);
% dEHeI   = reshape(Data.(,35),m,n);
% nuHeI   = reshape(Data.(,36),m,n);
nHeII   = reshape(Data.nHeII,m,n);
EHeII   = reshape(Data.EHeII,m,n); THeII=EHeII./(3/2*nHeII);
% dnHeII  = reshape(Data.(,39),m,n);
% dEHeII  = reshape(Data.(,40),m,n);
% nuHeII  = reshape(Data.(,41),m,n);
nHeIII  = reshape(Data.nHeIII,m,n);
EHeIII  = reshape(Data.EHeIII,m,n); THeIII=EHeIII./(3/2*nHeIII);
% dnHeIII = reshape(Data.(,44),m,n);
% dEHeIII  = reshape(Data.(,45),m,n);
% nuHeIII = reshape(Data.(,46),m,n);
% nCI     = reshape(Data.(,47),m,n);
% xnCI    = reshape(Data.(,48),m,n);
% nCII    = reshape(Data.(,49),m,n);
% xnCII   = reshape(Data.(,50),m,n);
% nCIII   = reshape(Data.(,51),m,n);
% xnCIII  = reshape(Data.(,52),m,n);
% nCIV    = reshape(Data.(,53),m,n);
% xnCIV   = reshape(Data.(,54),m,n);
% nCV     = reshape(Data.(,55),m,n);
% xnCV    = reshape(Data.(,56),m,n);
PRFe    = reshape(Data.PRFe,m,n);
PRFHi   = reshape(Data.PRFHi,m,n);
PRFH2i  = reshape(Data.PRFH2i,m,n);
PRFH3i  = reshape(Data.PRFH3i,m,n);
PRFHeII = reshape(Data.PRFHeII,m,n);
PRFHeIII= reshape(Data.PRFHeIII,m,n);
% dt      = reshape(Data.(:,63),m,n);
% dt = dt(1,:);

mycolor=jet(length(Rpos));
mycolor=lines(length(Rpos));
nel=0;
for i=1:length(Rpos)-1
    nel=nel+0.5*(ne(i+1,:)+ne(i,:))*(Rpos(i+1)^2-Rpos(i)^2)/(Rpos(end)^2-Rpos(1)^2);%
end

myidplot=round(length(Rpos)/2);

B=B0./Rpos*R;

%% plot data
if (video)
h=figure(20)
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

%     wce = interp(wce,10);
%     wcHi = interp(wcHi,10);
%     wcDi = interp(wcDi,10);
%     wcH2i = interp(wcH2i,10);
%     wcHDi = interp(wcHDi,10);
%     wcD2i = interp(wcD2i,10);
%     wcH3i = interp(wcH3i,10);
%     wcHeII = interp(wcHeII,10);
%     wcHeIII = interp(wcHeIII,10);
    Rpos_ = Rpos;%interp(Rpos,10);

    filename = [file(end-22:end-4)];

    for ind=[ 1:every:length(time)  length(time) ]
        wpe=5.64e4*sqrt(ne(:,ind));
        wpHi=1.32e3*Z*sqrt(nHi(:,ind)*HtoHD/1);
        wpDi=1.32e3*Z*sqrt(nHi(:,ind)*(1-HtoHD)/2);
        wpH2i=1.32e3*Z*sqrt(nH2i(:,ind)*HtoHD*HtoHD/2);
        wpHDi=1.32e3*Z*sqrt(nH2i(:,ind)*2*(1-HtoHD)*HtoHD/3);
        wpD2i=1.32e3*Z*sqrt(nH2i(:,ind)*(1-HtoHD)*(1-HtoHD)/4);
        wpH3i=1.32e3*Z*sqrt(nH3i(:,ind)/3);
        wpHeII=1.32e3*Z*sqrt(nHeII(:,ind)/4);
        wpHeIII=1.32e3*2*sqrt(nHeIII(:,ind)/4);

%         wpe = interp(wpe,10);
%         wpHi = interp(wpHi,10);
%         wpDi = interp(wpDi,10);
%         wpH2i = interp(wpH2i,10);
%         wpHDi = interp(wpHDi,10);
%         wpD2i = interp(wpD2i,10);
%         wpH3i = interp(wpH3i,10);
%         wpHeII = interp(wpHeII,10);
%         wpHeIII = interp(wpHeIII,10);

        wlh=(1 - ((wpe.^2)./((2*pi*f)^2-wce.^2)) ...
            - ((wpHi.^2)./((2*pi*f)^2-wcHi.^2)) ...
            - ((wpDi.^2)./((2*pi*f)^2-wcDi.^2)) ...
            - ((wpH2i.^2)./((2*pi*f)^2-wcH2i.^2)) ...
            - ((wpHDi.^2)./((2*pi*f)^2-wcHDi.^2)) ...
            - ((wpD2i.^2)./((2*pi*f)^2-wcD2i.^2)));% ...
            %- ((wpH3i.^2)./((2*pi*f)^2-wcH3i.^2)) ...
            %- ((wpHeII.^2)./((2*pi*f)^2-wcHeII.^2))  ...
            %- ((wpHeIII.^2)./((2*pi*f)^2-wcHeIII.^2)) );
        locs = (find(diff(sign(wlh))==2 | diff(sign(wlh))==-2));

        ax1=subplot(3,1,1)
        y=ne(:,ind);
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
        title(['Tomator-1D/WEST 48MHz HtoHD=' num2str(HtoHD) ',  ne(t=' num2str(time(ind)*1e3,'%10.2f') 'ms)'])
        legend('ne(r)','Antenna','LHR','Limiter','Location','NorthWest')
        ylabel('ne [1/cm3]')
%         xlabel('R pos [cm]')
        ylim([0 15e10])
        xlim([rmin rmax])
        set(gca,'FontSize',16) % right y-axis
        hold off

        ax2=subplot(3,1,2)
        y=Te(:,ind);
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
%         title(['Tomator-1D/KIPT-RF TOMAS H2 25MHz Te(t=' num2str(time(ind)*1e3,'%10.2f') 'ms)'])
        legend('Te(r)','Antenna','LHR','Limiter','Location','NorthWest')
        ylabel('Te [eV]')
%         xlabel('R pos [cm]')
        ylim([0 6])
        xlim([rmin rmax])
        set(gca,'FontSize',16)
        hold off

        ax2=subplot(3,1,3)
        plot(Rpos,wlh,'LineWidth',2); grid on; hold on
        ylim([-1 1])

%         y=([PRFe(:,ind)  PRFHi(:,ind)  PRFH2i(:,ind)  PRFH3i(:,ind)  PRFHeII(:,ind)  PRFHeIII(:,ind)])*1.6e-19;
%         plot(Rpos,y,'LineWidth',2); grid on; hold on
%         xline(radant,'.-','Color',mycolor(3,:),'LineWidth',4);
%         if (length(locs)>0)
%             xline(Rpos_(locs(1)),'--','Color',mycolor(2,:),'LineWidth',2);
%         end
% %         xline(Rpos_(locs(1)),'--','Color',mycolor(2,:),'LineWidth',2);
%         xline(R-lHFS*a,'k-');
%         xline(R+lLFS*a,'k-');
%         if (length(locs)>1)
%             for j=2:length(locs)
%                 xline(Rpos_(locs(j)),'--','Color',mycolor(2,:),'LineWidth',2);
%             end
%         end
% %         title('power deposition profiles(t_{end})')
%         legend('e','Hi','H2i','H3i','HeII','HeIII','Location','NorthWest')
%         ylabel('power [eV/s/cm3]')
%         ylim([0 0.04])
%         xlim([rmin rmax])
%         xlabel('R pos [cm]')
%         set(gca,'FontSize',16)
        hold off

        drawnow
        % Capture the plot as an image
        frame = getframe(h);
        im = frame2im(frame);
        if savevid
            imwrite(im,[filename '_' num2str(ind)],'png');
        end

    end
    linkaxes([ax1 ax2],'x')
end