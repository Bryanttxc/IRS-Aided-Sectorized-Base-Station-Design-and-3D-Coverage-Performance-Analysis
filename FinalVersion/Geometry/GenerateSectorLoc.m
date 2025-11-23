function [SectorLocation,LabelLocation,NumSectors] = GenerateSectorLoc()

InterSiteDist = 500; % m
r = InterSiteDist/3; % radius of the sector eNodeB cell
NumSiteRings = 1; % 0 ring->1 BS site, 1 ring->7 BS sites, 2 ring->19 BS sites, 3 ring->37 BS sites
rBorder = NumSiteRings*3*r+2*r; % radius of considered area

%%sector location and boresight
if NumSiteRings==0
    Sites=[0,0];%x,y coordinate of all Sites
elseif NumSiteRings==1
    Sites=zeros(7,2);
    ang=0:60:300;
    Sites(2:7,1)=3*r.*cosd(ang);
    Sites(2:7,2)=3*r.*sind(ang);
elseif NumSiteRings==2
    Sites=zeros(19,2);
    ang=0:60:300;
    Sites(2:7,1)=3*r.*cosd(ang);
    Sites(2:7,2)=3*r.*sind(ang);
    Sites(8:2:18,1)=6*r.*cosd(ang);
    Sites(8:2:18,2)=6*r.*sind(ang);  
    ang=30:60:330;
    Sites(9:2:19,1)=3*sqrt(3)*r.*cosd(ang);
    Sites(9:2:19,2)=3*sqrt(3)*r.*sind(ang);      
elseif NumSiteRings==3
    ;
end

Sites;
NumSites=size(Sites,1);

figure;
plot(Sites(:,1),Sites(:,2),'.')
hold on;
LabelLocation = zeros(3*NumSites,2);
linewidth = 1;
for p=1:NumSites
    %%%plot parallel lines: 60 degree lines
    plot([Sites(p,1) Sites(p,1)+r*cosd(60)],[Sites(p,2) Sites(p,2)+r*sind(60)],'-b','LineWidth',linewidth)
    plot([Sites(p,1) Sites(p,1)+r*cosd(60)]-1.5*r,[Sites(p,2) Sites(p,2)+r*sind(60)]+sqrt(3)/2*r,'-b','LineWidth',linewidth)
    plot([Sites(p,1) Sites(p,1)+r*cosd(60)]+1.5*r,[Sites(p,2) Sites(p,2)+r*sind(60)]-sqrt(3)/2*r,'-b','LineWidth',linewidth)
    plot([Sites(p,1) Sites(p,1)+r*cosd(60)]-1.5*r,[Sites(p,2) Sites(p,2)+r*sind(60)]-sqrt(3)/2*r,'-b','LineWidth',linewidth)
    plot([Sites(p,1) Sites(p,1)+r*cosd(60)],[Sites(p,2) Sites(p,2)+r*sind(60)]-sqrt(3)*r,'-b','LineWidth',linewidth)

    %%%plot parallel lines: 180 degree lines
    plot([Sites(p,1) Sites(p,1)+r*cosd(180)],[Sites(p,2) Sites(p,2)+r*sind(180)],'-b','LineWidth',linewidth)
    plot([Sites(p,1) Sites(p,1)+r*cosd(180)],[Sites(p,2) Sites(p,2)+r*sind(180)]+sqrt(3)*r,'-b','LineWidth',linewidth)
    plot([Sites(p,1) Sites(p,1)+r*cosd(180)],[Sites(p,2) Sites(p,2)+r*sind(180)]-sqrt(3)*r,'-b','LineWidth',linewidth)
    plot([Sites(p,1) Sites(p,1)+r*cosd(180)]+1.5*r,[Sites(p,2) Sites(p,2)+r*sind(180)]+sqrt(3)/2*r,'-b','LineWidth',linewidth)
    plot([Sites(p,1) Sites(p,1)+r*cosd(180)]+1.5*r,[Sites(p,2) Sites(p,2)+r*sind(180)]-sqrt(3)/2*r,'-b','LineWidth',linewidth)
    
    %%%plot parallel lines: 300 degree lines  
    plot([Sites(p,1) Sites(p,1)+r*cosd(300)],[Sites(p,2) Sites(p,2)+r*sind(300)],'-b','LineWidth',3)
    plot([Sites(p,1) Sites(p,1)+r*cosd(300)]-1.5*r,[Sites(p,2) Sites(p,2)+r*sind(300)]-sqrt(3)/2*r,'-b','LineWidth',linewidth)
    plot([Sites(p,1) Sites(p,1)+r*cosd(300)]+1.5*r,[Sites(p,2) Sites(p,2)+r*sind(300)]+sqrt(3)/2*r,'-b','LineWidth',linewidth)
    plot([Sites(p,1) Sites(p,1)+r*cosd(300)]-1.5*r,[Sites(p,2) Sites(p,2)+r*sind(300)]+sqrt(3)/2*r,'-b','LineWidth',linewidth)
    plot([Sites(p,1) Sites(p,1)+r*cosd(300)],[Sites(p,2) Sites(p,2)+r*sind(300)]+sqrt(3)*r,'-b','LineWidth',linewidth)
    
    %%%plot the sector arrow
    arrow3(Sites(p,:),Sites(p,:)+[0.4*r 0],[],2.5,1.5);%s,W,H  see arrow3 for their meaning
    arrow3(Sites(p,:),Sites(p,:)+0.4*r*[cosd(120) sind(120)],[],2.5,1.5);%s,W,H  see arrow3 for their meaning
    arrow3(Sites(p,:),Sites(p,:)+0.4*r*[cosd(240) sind(240)],[],2.5,1.5);%s,W,H  see arrow3 for their meaning

    %%%label the three sectors at each site
    LabelLocation(3*(p-1)+1,:) = [Sites(p,1)+r,Sites(p,2)];
    LabelLocation(3*(p-1)+2,:) = [Sites(p,1)-0.5*r,Sites(p,2)+sqrt(3)/2*r];
    LabelLocation(3*(p-1)+3,:) = [Sites(p,1)-0.5*r,Sites(p,2)-sqrt(3)/2*r];
    text(Sites(p,1)+0.75*r,Sites(p,2),sprintf('%d',3*(p-1)+1),'FontSize',28,'FontWeight','bold');
    text(Sites(p,1)-0.75*r,Sites(p,2)+0.75*r,sprintf('%d',3*(p-1)+2),'FontSize',28,'FontWeight','bold');
    text(Sites(p,1)-0.75*r,Sites(p,2)-0.75*r,sprintf('%d',3*(p-1)+3),'FontSize',28,'FontWeight', 'bold');
end
axis equal tight
xlim([-rBorder,rBorder]);
ylim([-rBorder,rBorder]);
xlabel('x (m)','Fontsize',30);
ylabel('y (m)','Fontsize',30);
set(gcf,'Color','w')
set(gca,'FontSize',15);
axis image
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'Visible','off') %// handy for printing the image
%axis equal

%%% sign: site %%%
hold on
plot(Sites(1,1)+1.4*r,Sites(1,2)+0.3*r,'k^','MarkerFaceColor','r','MarkerSize',20);
text(Sites(4,1)-2.5*r,Sites(4,2)+0.65*r,sprintf('站点'),'FontSize',28,'Color','r','FontWeight','bold')
text(Sites(4,1)-2.5*r,Sites(4,2)+0.60*r,sprintf('用户'),'FontSize',28,'Color','r','FontWeight','bold')
% for mainSignalPowerScalingLawinSingleCell.m
% hold on
% plot(60,-4.3,'b^','MarkerFaceColor','b','MarkerSize',10); % UEidx = 500

NumSectors = 3*NumSites;
SectorLocation = zeros(NumSectors,2);
for p = 1:NumSites
    SectorLocation(3*(p-1)+1:3*p,:) = repmat( Sites(p,:),[3 1]);
end
SectorLocation;
Boresight = zeros(NumSectors,1); %%degree
Boresight(2:3:end) = 120;
Boresight(3:3:end) = -120;
Boresight;

