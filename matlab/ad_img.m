%  This program draws the angular distribution and the image.
%  
%  
clf
%sets font size
set (gcf, 'DefaultTextFontSize',14)
set (gcf, 'DefaultAxesFontSize',14)
set(gcf,'PaperPosition',[1.125 2.52 6.25 5.96]);
%

colormap(gray);
format bank;
%Starting_time=clock

%-------------------------------
x=linspace(-1,1,128);
y=linspace(-4,4,128);
[X,Y]=meshgrid(x,y);
Z=2.4*exp(-X.^2/(0.04^2)-Y.^2/(0.15^2)); 
%------------------------------
e=1.6e-19;
eps=0.664e-19;                     %540nm...0.4149ev = 0.664e-19 J
V=2.35;  %(2.35 for 2.2v)
d=0.5;
D=1.35;
A=e*V*d/eps/D;
rmax=2*eps*D*sqrt(1+ A)/e/V;

%----------------------------
xscale = 45;
yscale =47;
for k=1:128
    for j=1:128
	xc=(k-64.5)/xscale;
	yc=(j-64.5)/yscale;
	x1 = abs(xc) - 0.5/xscale;
        x4 = abs(xc) + 0.5/xscale;
        y1 = abs(yc) - 0.5/xscale;
        y4 = abs(yc) + 0.5/xscale;
        r1=sqrt(x1^2+y1^2);
        r4=sqrt(x4^2+y4^2);
	 if r1>=rmax
	    Ip(k,j)=0;
	    projectionp(k,j) = 0;
	    projectionn(k,j) = 0;
	    mask(k,j)=0;
	 else
	    if r4<=rmax;
	       area = 1;
	    else
	   	r2=sqrt(x1^2+y4^2);
        	r3=sqrt(x4^2+y1^2);
        	   if r2>=rmax;
		      if r3>=rmax
		      	y2 = sqrt(rmax^2-x1^2);
		       	x3 = sqrt(rmax^2-y1^2);
                        area = 0.5*(y2-y1)*(x3-x1)*xscale*yscale;
		      	xc = sign(xc)*(2*x1/3 + x3/3);
			yc = sign(yc)*(2*y1/3 + y2/3);
		      else
		      	y2 = sqrt(rmax^2-x1^2);
                        y4 = sqrt(rmax^2-x4^2);
                        area = ((y2 + y4)/2 - y1)*yscale;
			yc = sign(yc)*((y2 + y4)/4 + y1/2);
		      end
		   else
                     if r3>=rmax;
		        x3 = sqrt(rmax^2-y1^2);
                    	x4 = sqrt(rmax^2-y4^2);
                    	area = ((x3 + x4)/2 - x1)*xscale;
                    	xc = sign(xc)*((x3 + x4)/4 + x1/2);
		     else
		        y6 = sqrt(rmax^2-x4^2);
                        x5 = sqrt(rmax^2-y4^2);
                        area = 1 - 0.5*(y4-y6)*(x4-x5)*xscale*yscale;
                        xc = (xc - sign(xc)*(2*x4/3 + x5/3)*(1-area))/area;
                        yc = (yc - sign(yc)*(2*y4/3 + y6/3)*(1-area))/area;
		     end
	          end
	    end
	    if xc>0
		phi = atan(yc/xc);
	    else
		phi = pi + atan(yc/xc);
	    end

  %--------------------------------------
       r = sqrt(xc^2 + yc^2);
	    B=e*V*r/eps/D;
	    C=A/B;
	    R=sqrt(1+A-(B/2)^2);
	    if B>sqrt(A)*2
	      	thetap = pi - asin(sqrt((1+A/2+R)/2/(C^2+1)));
	    else
		thetap= asin(sqrt((1+A/2+R)/2/(C^2+1)));
	    end
            thetan = pi - asin(sqrt((1+A/2-R)/2/(C^2+1)));
	    M=B^2/16/(C^2+1)/r^2/R;
	    Lp=C^2*(1+A/2+R)/2/(C^2+1)^2/r^2;
            Ln=C^2*(1+A/2-R)/2/(C^2+1)^2/r^2;
	    costhap=sin(thetap)*cos(phi);
            costhan=sin(thetan)*cos(phi);
	    if thetap<pi/2
	        phiap = pi-atan(tan(thetap)*sin(phi));
	    else
	        phiap = -atan(tan(thetap)*sin(phi));
	    end
	    phian = -atan(tan(thetan)*sin(phi));
	    projectionp(k,j)=abs((Lp-M)/cos(thetap))*area;
	    projectionn(k,j)=abs((Ln+M)/cos(thetan))*area;
	    p2=legendre(2,costhap);
            Y20p(k,j) = sqrt(5/4/pi)*p2(1);
            Y21p(k,j) = sqrt(5/4/pi/6)*p2(2)*exp(i*phiap);  %Signs differ in Matlab and Arfkin
            Y22p(k,j) = sqrt(5/4/pi/24)*p2(3)*exp(2*i*phiap);
%            Y2m1p(k,j) = -1*conj(Y21);
%            Y2m2p(k,j) = conj(Y22);

 	    p2=legendre(2,costhan);
            Y20n(k,j) = sqrt(5/4/pi)*p2(1);
            Y21n(k,j) = sqrt(5/4/pi/6)*p2(2)*exp(i*phian);
            Y22n(k,j) = sqrt(5/4/pi/24)*p2(3)*exp(2*i*phian);
%            Y2m1n(k,j) = -1*conj(Y21);
%            Y2m2n(k,j) = conj(Y22);
     end
  end
end

            Y00 = 1/sqrt(4*pi);


polariz = [-1   0.41*i];
label=['abcd'];
AOB = 0.38;
COB = -1.05;
del = 5.08;

for ll = 1:2
 pol = polariz(ll);
%
phiap=(0.:0.1*pi:2*pi);
thetap1=(0.00:.03*pi:pi);
[THETA,PHI] = meshgrid(thetap1,phiap);
costhap=cos(THETA);
          p2b=legendre(2,costhap);
            Y20pb(:,:) = sqrt(5/4/pi).*p2b(1,:,:);
            pp2b(:,:)=p2b(2,:,:);
            Y21pb(:,:) = sqrt(5/4/pi/6)*pp2b(:,:).*exp(i.*PHI);  %Signs differ in %Matlab and Arfkin
            pp3b(:,:)=p2b(3,:,:);
            Y22pb = sqrt(5/4/pi/24)*pp3b(:,:).*exp(2*i*PHI);

   shpuub=(abs(Y00*exp(i*del)*AOB*sqrt(10/3)*(1+pol^2) + ...
     (Y20pb*(-2+pol^2)*sqrt(2/3) + 4*pol*real(Y21pb) - ...
       2*pol^2*real(Y22pb)) + ... 
     COB/3*(-pol*2*i*imag(Y21pb) + pol^2 * 2*i*imag(Y22pb)))).^2;
   shpudb = (COB/3)^2 * (abs(2*Y21pb + pol*(sqrt(6)*Y20pb - ...
      2*Y22pb)-...
       2* pol^2 * real(Y21pb))).^2;
   shpddb=(abs(Y00*exp(i*del)*AOB*sqrt(10/3)*(1+pol^2) + ...
     (Y20pb*(-2+pol^2)*sqrt(2/3) + 4*pol*real(Y21pb) - ...
        2*pol^2*real(Y22pb)) + ... 
      COB/3*(pol*2*i*imag(Y21pb) - pol^2 * 2*i*imag(Y22pb)))).^2;
   shpdub = (COB/3)^2 * (abs(-2*conj(Y21pb) + ...
      pol*(-sqrt(6)*Y20pb + ...
        2*conj(Y22pb)) +    2*pol^2 * real(Y21pb))).^2;
   shpb = (shpuub + shpudb + shpddb + shpdub)/(1+(abs(pol))^2).^2;


subplot(2,2,ll)

[XX,YY,ZZ]=sph2cart(PHI,pi/2-THETA,shpb);
mesh(-ZZ,YY,XX)
%set(gcf,'Color',[0 0 0])
view(0,0)
hold on
axis([-2 2 -2 2 -2 2])
text(-1.8,0,1.5,['(',label(ll),')'])
set(gca,'Visible','off')

 for k = 1:128
 for j = 1:128
  if projectionp(k,j) ~= 0
   shpuu=(abs(Y00*exp(i*del)*AOB*sqrt(10/3)*(1+pol^2) + ...
     (Y20p(k,j)*(-2+pol^2)*sqrt(2/3) + 4*pol*real(Y21p(k,j)) - ...
       2*pol^2*real(Y22p(k,j))) + ... 
     COB/3*(-pol*2*i*imag(Y21p(k,j)) + pol^2 * 2*i*imag(Y22p(k,j)))))^2;
   shpud = (COB/3)^2 * (abs(2*Y21p(k,j) + pol*(sqrt(6)*Y20p(k,j) - ...
      2*Y22p(k,j))-...
       2* pol^2 * real(Y21p(k,j))))^2;
   shpdd=(abs(Y00*exp(i*del)*AOB*sqrt(10/3)*(1+pol^2) + ...
     (Y20p(k,j)*(-2+pol^2)*sqrt(2/3) + 4*pol*real(Y21p(k,j)) - ...
        2*pol^2*real(Y22p(k,j))) + ... 
      COB/3*(pol*2*i*imag(Y21p(k,j)) - pol^2 * 2*i*imag(Y22p(k,j)))))^2;
   shpdu = (COB/3)^2 * (abs(-2*conj(Y21p(k,j)) + ...
      pol*(-sqrt(6)*Y20p(k,j) + ...
        2*conj(Y22p(k,j))) +    2*pol^2 * real(Y21p(k,j))))^2;
   shp = (shpuu + shpud + shpdd + shpdu)/(1+(abs(pol))^2)^2;
     
   shpuu=(abs(Y00*exp(i*del)*AOB*sqrt(10/3)*(1+pol^2) + ...
     (Y20n(k,j)*(-2+pol^2)*sqrt(2/3) + 4*pol*real(Y21n(k,j)) - ...
        2*pol^2*real(Y22n(k,j))) + ... 
    COB/3*(-pol*2*i*imag(Y21n(k,j)) + pol^2 * 2*i*imag(Y22n(k,j)))))^2;
   shpud = (COB/3)^2 * (abs(2*Y21n(k,j) + pol*(sqrt(6)*Y20n(k,j) - ...
       2*Y22n(k,j))-...
         2* pol^2 * real(Y21n(k,j))))^2;
   shpdd=(abs(Y00*exp(i*del)*AOB*sqrt(10/3)*(1+pol^2) + ...
          (Y20n(k,j)*(-2+pol^2)*sqrt(2/3) + 4*pol*real(Y21n(k,j)) - ...
               2*pol^2*real(Y22n(k,j))) + ... 
          COB/3*(pol*2*i*imag(Y21n(k,j)) - pol^2 * 2*i*imag(Y22n(k,j)))))^2;
   shpdu = (COB/3)^2 * (abs(-2*conj(Y21n(k,j)) + ...
      pol*(-sqrt(6)*Y20n(k,j) + ...
          2*conj(Y22n(k,j))) +    2*pol^2 * real(Y21n(k,j))))^2;
   shn = (shpuu + shpud + shpdd + shpdu)/(1+(abs(pol))^2)^2;
      
   Ip(k,j)=shp*projectionp(k,j)+ shn*projectionn(k,j);
   end;
  end; 
end;

Ic=conv2(Ip,Z,'same');
subplot(2,2,ll+2)
pcolor(Ic(1:2:128,:)')
shading flat;
hold on
text(5,115,['(',label(ll+2),')'],'Color',[1 1 1])

end


