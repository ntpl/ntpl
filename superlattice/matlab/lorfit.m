function widthhm = lorfit(w,b)
	PT_PERC=1E-3;
   	 PEAK_PERC = 1.5;
    	gamma_guess = 0.5;
	[Ipeak,Jpeak]=max(b);
    	INV_PERC=0.1;

	[I,J] = find(b(1:Jpeak)<PT_PERC*b(Jpeak));
    %I
    if length(I)<1
        wleft = 1;
    else
        wleft = I(length(I));
    end
	[I,J] = find(b(Jpeak:length(w))>PT_PERC*b(Jpeak));
	wright = Jpeak + I(length(I));
	if wright > length(b)
	 	wright=length(b);
	end
    
    	weights = ones(length(w(wleft:wright)),1);
    	weights(1:5) = INV_PERC/PT_PERC;
    	weights(length(weights)-5:length(weights)) = INV_PERC/PT_PERC;
    
 	lor_func = @(c,w)weights.*(c(1))./(1 + ((w - c(3))./ c(2) ).^2 );
	options = optimset('MaxIter',10000,'MaxFunEvals',10000,'ScaleProblem','Jacobian','TolX',1e-6,'TolFun',1e-6);
	c0 = [ PEAK_PERC*Ipeak, gamma_guess, w(Jpeak)];
	lb(1:length(c0)) = 0.0; 
        %ub(1:3:length(c0)) = PEAK_PERC*Ipeak; 
        ub(1:3:length(c0)) = max(b)*10; 
        ub(2:3:length(c0)) = 1000; 
        ub(3:3:length(c0)) = w(length(w));
	[c_fit,resnorm,residual] = lsqcurvefit(lor_func,c0,w(wleft:wright),b(wleft:wright).*weights,lb,ub,options);
	resnorm
    	widthhm=c_fit(2)
	%dummy=lor_func;
    %semilogy(w(wleft:wright),lor_func(c0,w(wleft:wright))./weights)
	%hold on
	%semilogy(w(wleft:wright),lor_func(c_fit,w(wleft:wright))./weights)
	%hold on
	%semilogy(w,b,'ro')
	%hold off
end
