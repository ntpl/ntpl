function plotwvk

NMD.x0.LJ.eps = 1.67E-21;              
NMD.x0.LJ.sigma = 3.4E-10;
NMD.x0.LJ.mass = 6.6326E-26;
NMD.x0.LJ.tau = sqrt((NMD.x0.LJ.mass*(NMD.x0.LJ.sigma)^2)/NMD.x0.LJ.eps);
kb = 1.3806E-23; 

p4_freq='4p_freq.dat'
p8_freq='8p_freq.dat'
p12_freq='12p_freq.dat'
p16_freq='16p_freq.dat'

p4_vel='4p_vel.dat'
p8_vel='8p_vel.dat'
p12_vel='12p_vel.dat'
p16_vel='16p_vel.dat'

p4_life='4p_life.dat'
p8_life='8p_life.dat'
p12_life='12p_life.dat'
p16_life='16p_life.dat'

p4_x0='4p_x0.dat'
p8_x0='8p_x0.dat'
p12_x0='12p_x0.dat'
p16_x0='16p_x0.dat'

p6_freq='6p_freq.dat'
p6_vel='6p_vel.dat'
p6_life='6p_life.dat'
p6_x0='6p_x0.dat'

p10_freq='10p_freq.dat'
p10_vel='10p_vel.dat'
p10_life='10p_life.dat'
p10_x0='10p_x0.dat'

[p4f,p4_kw] =omegabin(p4_freq,p4_vel,p4_life,p4_x0);
[p6f,p6_kw] =omegabin(p6_freq,p6_vel,p6_life,p6_x0);
[p8f,p8_kw] =omegabin(p8_freq,p8_vel,p8_life,p8_x0);
[p10f,p10_kw] =omegabin(p10_freq,p10_vel,p10_life,p10_x0);
[p12f,p12_kw] =omegabin(p12_freq,p12_vel,p12_life,p12_x0);
[p16f,p16_kw] =omegabin(p16_freq,p16_vel,p16_life,p16_x0);

semilogx(p4f,p4_kw,'bo',p6f,p6_kw,'co',p8f,p8_kw,'go',p10f,p10_kw,'mo',p12f,p12_kw,'ro',p16f,p16_kw,'ko')
legend('4p','6p','8p','10p','12p','16p');
xlabel('Mean Free Path [m]')
ylabel('Thermal Conductivity [W/m K]')

yL = get(gca,'YLim');
hold on
line([4*2*0.78161*NMD.x0.LJ.sigma 4*2*0.78161*NMD.x0.LJ.sigma],yL,'Color','b');
line([6*2*0.78161*NMD.x0.LJ.sigma 6*2*0.78161*NMD.x0.LJ.sigma],yL,'Color','c')
line([8*2*0.78161*NMD.x0.LJ.sigma 8*2*0.78161*NMD.x0.LJ.sigma],yL,'Color','g')
line([10*2*0.78161*NMD.x0.LJ.sigma 10*2*0.78161*NMD.x0.LJ.sigma],yL,'Color','m')
line([12*2*0.78161*NMD.x0.LJ.sigma 12*2*0.78161*NMD.x0.LJ.sigma],yL,'Color','r')
line([16*2*0.78161*NMD.x0.LJ.sigma 16*2*0.78161*NMD.x0.LJ.sigma],yL,'Color','k')

end

function [fkw,kw]=omegabin(str_freq,str_vel,str_life,str_vol)

NMD.x0.LJ.eps = 1.67E-21;              
NMD.x0.LJ.sigma = 3.4E-10;
NMD.x0.LJ.mass = 6.6326E-26;
NMD.x0.LJ.tau = sqrt((NMD.x0.LJ.mass*(NMD.x0.LJ.sigma)^2)/NMD.x0.LJ.eps);
kb = 1.3806E-23; 

ff=reshape(load(str_freq)',[],1);
vel=load(str_vel)*(NMD.x0.LJ.sigma/NMD.x0.LJ.tau);
lifetime=load(str_life);
x0=load(str_vol);

NUM_ATOMS=x0(1,1);
L(1) = x0(1,3); L(2) = x0(1,4); L(3) = x0(1,5); 
VOLUME = (L(1)*L(2)*L(3)*NMD.x0.LJ.sigma^3);

velx=reshape(vel(:,1),size(lifetime,2),size(lifetime,1))';
vely=reshape(vel(:,2),size(lifetime,2),size(lifetime,1))';
velz=reshape(vel(:,3),size(lifetime,2),size(lifetime,1))';

kappax = sum(sum((kb/VOLUME).*lifetime.*((velx).^2)))
kappay = sum(sum((kb/VOLUME).*lifetime.*(vely.^2)))
kappaz = sum(sum((kb/VOLUME).*lifetime.*(velz.^2)))
khs = 3/2*(pi/6)^(1/3)*kb*(NUM_ATOMS/VOLUME)^(2/3)*(0.8*max(reshape(velx.',[],1)))

ll=reshape(lifetime.',[],1);
vx=reshape(velx.',[],1);
vy=reshape(vely.',[],1);
vz=reshape(velz.',[],1);


%m(:,1)=ff;
m(:,1)=sqrt(vel(:,1).^2+vel(:,2).^2+vel(:,3).^2).*ll;
m(:,2)=(kb/VOLUME).*ll.*(vx.^2);
%m(:,2)=((kb/VOLUME).*ll.*(vy.^2)+(kb/VOLUME).*ll.*(vz.^2))/2;


m=sortrows(m,1);
indices = find(m(:,1)==0);
m(indices,:) = [];

kw=zeros(200,1);
dw=floor(length(m(:,1))/length(kw));
fdw=max(m(:,1))/length(kw);
fkw=fdw*(1:1:length(kw));

for j=2:1:length(kw)
    [I]=logical(fdw*(j-1)<m(:,1) & m(:,1)<fdw*j );
    kw(j)=kw(j-1)+sum(m(I,2));
end

%semilogx(fdw*(1:1:length(kw)),kw,'o')
end



