/*
Name:	lmetapreg.do
Creator:	Victoria N Nyaga
Date: 19th October 2018
Purpose:	Facilitate maintanance of mata functions in lmetaprop.mlib.
			These functions are used in metapreg.
*/

cd "C:\ado\plus\l"
mata:mata mlib create lmetapreg, replace

version 14.0 
mata:
	mata clear

//************************************************************************ Added 7.12.2020
	void hetgrp(real rowvector ES, real rowvector se) {
			real scalar mu, chisq, g, diff
			g = length(se)
			ones = J(1, g, 1)
			
			mu = (sum(ES :/ (se :* se)))/ sum(ones :/ (se :* se))
			MU = J(1, g, mu)
			diff = (ES - MU)
			chisq = sum((diff :* diff) :*(ones :/ (se :* se)))
			st_local("chisq", chisq)
	}

//************************************************************************
	void bayesRR(real rowvector x, real scalar nsim, real scalar seed) {

		rseed(seed)
		real scalar a1, b1, a2, b2, prop1, prop2
		RR = J(nsim, 1, .)
		
		a1 = x[1]
		b1 = x[2]
		a2 = x[3]
		b2 = x[4]
		
		for (i=1; i<=nsim; i++) {
			prop1 = rbeta(1, 1, a1, b1)
			prop2 = rbeta(1, 1, a2, b2)		
			RR[i] = prop1/prop2
		}

		st_matrix("RR", RR)
	}
//************************************************************************	
	void cml_ci(real rowvector x, real scalar alpha) {
	
	real scalar a, b, c, d, n, p1, p0, RR, inits_l, inits_u

	ci = J(1, 2, .)

	a = x[1]
	b = x[2]
	c = x[3]
	d = x[4]
	
	n = a + b + c + d
	
	p1 = (a + b)/n
	p0 = (a + c)/n
	
	RR = p1/p0
	
	inits_l = RR*exp(-1*invnormal(1 - alpha/2)*sqrt((b + c)/((a + b)*(a + c))))*0.9
	inits_u = RR*exp(invnormal(1 - alpha/2)*sqrt((b + c)/((a + b)*(a + c))))*1.1	
		
	S = solvenl_init()
	solvenl_init_evaluator(S, &cmlfun())
	solvenl_init_type(S, "zero")
	solvenl_init_technique(S, "newton")
	solvenl_init_numeq(S, 1)
	solvenl_init_narguments(S, 1)	
	
	//solvenl_init_argument(S, 1, (a, b, c, n, alpha, 1))
	solvenl_init_argument(S, 1, (a, b, c, n, alpha))
	solvenl_init_startingvals(S, inits_l)
	solvenl_init_iter_log(S, "off")
	lower = solvenl_solve(S)

	//solvenl_init_argument(S, 1, (a, b, c, n, alpha, -1))
	solvenl_init_startingvals(S, inits_u)
	solvenl_init_iter_log(S, "off")
	upper = solvenl_solve(S)
	
	ci = (lower, upper)
	st_matrix("ci", ci)
}
//************************************************************************
//score
//	void fun(real scalar RR, real scalar y, real rowvector x){
//		a = x[1]
//		b = x[2]
//		c = x[3]
//		n = x[4]
//		alpha = x[5]
//		sign = x[6]
//		
//		real scalar p10_tilde
//		real scalar p01_tilde
//		
//		z_alpha = invnormal(alpha/2)
//		
//		A = n*(1 + RR)
//		B = (c + a)*RR^2 - (a + b + 2*c)
//		C = c*(1 - RR)*(a + b + c)/n
//		p01_tilde = (-B + sqrt(B^2 - 4*A*C))/(2*A)
//		p10_tilde = ((a + b + c)/n - p01_tilde)/RR
//		y = ((b + a - (a + c)*RR)/(RR*(2*p01_tilde + p10_tilde*(RR - 1))))/sqrt((n)/(RR*(2*p01_tilde + p10_tilde*(RR - 1)))) + sign*z_alpha
//	}
//************************************************************************
	void cmlfun(real scalar RR, real scalar y, real rowvector df){
		n11 = df[1] 
		n10 = df[2]
		n01 = df[3]
		n00 = df[4]
		alpha = df[5]

		n = n11 + n01 + n10 + n00
		
		p1_hat = (n11 + n10)/n
		p0_hat = (n11 + n01)/n
		p01_hat = n01/n
		p10_hat = n10/n
		p00_hat = n00/n       

		p10_tilde =(-p1_hat + (RR^2)*(p0_hat + 2*p10_hat) + sqrt((p1_hat - (RR^2)*p0_hat)^2 + 4*(RR^2)*p10_hat*p01_hat))/(2*RR*(RR + 1))
		p01_tilde = RR*p10_tilde - (RR - 1)*(1 - p00_hat)
		y = n*((p1_hat - RR*p0_hat)^2)/(RR*(p10_tilde + p01_tilde)) - invnormal(alpha/2)^2
	}
//************************************************************************	
	void koopmancifun(real scalar RR, real scalar y, real rowvector x){
		x1 = x[1]
		n1 = x[2]
		x2 = x[3]
		n2 = x[4]
		level = x[5]
	
		real scalar p1_tilde
		
		ki = invchi2(1, 1 - level)		
		
		p1_tilde = (RR*(n1 + x2) + x1 + n2 - sqrt((RR*(n1 + x2) + x1 + n2)^2 - 4*RR*(n1 + n2)*(x1 + x2)))/(2*(n1 + n2))
		y = (((x1 - n1*p1_tilde)^2)/(n1*p1_tilde*(1 - p1_tilde)))*(1 + (n1*(RR - p1_tilde))/(n2*(1 - p1_tilde))) -  ki
	}

//************************************************************************
	void bootpi_matchre(real rowvector x, real scalar breps, real scalar seed){		
		
		rseed(seed)
		real scalar b0, b1, s_b0, v_b0_plus_b1, tausq, w_tausq
		predodds = J(breps, 2, .)
		predRR = J(breps, 1, .)
		
		b0 = x[1]
		b1 = x[2]
		v_b0 = x[3]
		v_b0_plus_b1 = x[4]
		tausq = x[5]
		w_tausq = x[6]		
		
		for (i=1; i<=breps; i++) {
			predodds[i, 1] = rnormal(1, 1, b0, sqrt(w_tausq + tausq + v_b0))
			predodds[i, 2] = rnormal(1, 1, b0 + b1, sqrt(w_tausq + tausq + v_b0_plus_b1))
			predRR[i, 1] = invlogit(predodds[i, 2])/invlogit(predodds[i, 1])						
		}
		st_matrix("predRR", predRR)
	}

//************************************************************************
	void bootpi_matchfe(real rowvector x, real scalar breps, real scalar seed){		
		
		rseed(seed)
		real scalar b0, b1, v_b0, v_b0_plus_b1, tausq
		predodds = J(breps, 2, .)
		predRR = J(breps, 1, .)
		
		b0 = x[1]
		b1 = x[2]
		v_b0 = x[3]
		v_b0_plus_b1 = x[4]
		tausq = x[5]		
		
		for (i=1; i<=breps; i++) {
			predodds[i, 1] = rnormal(1, 1, b0, sqrt(tausq + v_b0))
			predodds[i, 2] = rnormal(1, 1, b0 + b1, sqrt(tausq + v_b0_plus_b1))
			predRR[i, 1] = invlogit(predodds[i, 2])/invlogit(predodds[i, 1])						
		}
		st_matrix("predRR", predRR)
	}

//************************************************************************
//Added on 17 Jan 2019
void koopman_ci(real rowvector v, real scalar alpha) {
	real scalar zstar, x, m, y, n, CIL, CIU, nrat, varhat
	
	zstar = invnormal(1 - alpha/2)
	x = v[1] 
	m = v[2] 
	y = v[3] 
	n = v[4] 
	
	if (x == 0 & y == 0) {
		CIL = 0
		CIU = .
	}
	else {
	    a1 = n * (n * (n + m) * x + m * (n + x) * (zstar^2))
        a2 = -n * (n * m * (y + x) + 2 * (n + m) * y * x + m * (n + y + 2 * x) * (zstar^2))
        a3 = 2 * n * m * y * (y + x) + (n + m) * (y^2) * x + n * m * (y + x) * (zstar^2)
        a4 = -m * (y^2) * (y + x)
        b1 = a2/a1
        b2 = a3/a1
        b3 = a4/a1
        c1 = b2 - (b1^2)/3
        c2 = b3 - b1 * b2/3 + 2 * (b1^3)/27
        ceta = (acos(sqrt(27) * c2/(2 * c1 * sqrt(-c1))))
        t1 = (-2 * sqrt(-c1/3) * cos(pi()/3 - ceta/3))
        t2 = (-2 * sqrt(-c1/3) * cos(pi()/3 + ceta/3))
        t3 = (2 * sqrt(-c1/3) * cos(ceta/3))
        p01 = t1 - b1/3
        p02 = t2 - b1/3
        p03 = t3 - b1/3
        p0sum = p01 + p02 + p03
        p0up = min((p01, p02, p03))
        p0low = p0sum - p0up - max((p01, p02, p03))
        
        rat = (x/m)/(y/n)
        nrat = (x/m)/(y/n)
        varhat = (1/x) - (1/m) + (1/y) - (1/n)
        
        if ((x == 0) & (y != 0)) {
            nrat = ((x + 0.5)/m)/(y/n)
            varhat = (1/(x + 0.5)) - (1/m) + (1/y) - (1/n)
        }
        if ((y == 0) & (x != 0)) {
            nrat = (x/m)/((y + 0.5)/n)
            varhat = (1/x) - (1/m) + (1/(y + 0.5)) - (1/n)
        }
        if ((y == n) & (x == m)) {
            nrat = 1
            varhat = (1/(m - 0.5)) - (1/m) + 1/(n - 0.5) - (1/n)
        }
        La = nrat * exp(-1 * zstar * sqrt(varhat)) * 1/4
        Ha = nrat
	
	if ((x != 0) & (y == 0)) {
	  if (x == m) {
		CIL = (1 - (m - x) * (1 - p0low)/(y + m - (n + m) * p0low))/p0low
		CIU = .
	  }
	  else {
		S = solvenl_init()
		solvenl_init_evaluator(S, &U())
		solvenl_init_type(S, "zero")
		solvenl_init_technique(S, "newton")
		solvenl_init_numeq(S, 1)
		solvenl_init_narguments(S, 1)	
		solvenl_init_argument(S, 1, (x, m, y, n, alpha))
		solvenl_init_startingvals(S, La)
		solvenl_init_iter_log(S, "off")
		
		CIL = solvenl_solve(S)
		CIU = .
	  }
	}
	
	if ((x == 0) & (y != n)) {
		S = solvenl_init()
		solvenl_init_evaluator(S, &U())
		solvenl_init_type(S, "zero")
		solvenl_init_technique(S, "newton")
		solvenl_init_numeq(S, 1)
		solvenl_init_narguments(S, 1)	
		solvenl_init_argument(S, 1, (x, m, y, n, alpha))
		solvenl_init_startingvals(S, Ha)
		solvenl_init_iter_log(S, "off")
		
		CIU = solvenl_solve(S)
		CIL = 0
	}
	
	if (((x == m) | (y == n)) & (y != 0)) {
	  if ((x == m) & (y == n)) {
		CIL = m/(m + invchi2(1, 1 - alpha))
		CIU = (n + invchi2(1, 1 - alpha))/n	
	  }
	  if ((x == m) & (y != n)) {
		phat1 = x/m
		phat2 = y/n
		phihat = phat2/phat1
		phiu = 1.1 * phihat
		r = 0
		while (r >= -zstar) {
		  a = (m + n) * phiu
		  b = -((x + n) * phiu + y + m)
		  c = x + y
		  p1hat = (-b - sqrt(b^2 - 4 * a * c))/(2 * a)
		  p2hat = p1hat * phiu
		  q2hat = 1 - p2hat
		  var = (m * n * p2hat)/(n * (phiu - p2hat) + m * q2hat)
		  r = ((y - n * p2hat)/q2hat)/sqrt(var)
		  phiu1 = phiu
		  phiu = 1.0001 * phiu1
		}
		CIU = (1 - (m - x) * (1 - p0up)/(y + m - (n + m) * p0up))/p0up
		CIL = 1/phiu1
	  }
	  
	  if ((y == n) & (x != m)) {
		phat2 = y/n
		phat1 = x/m
		phihat = phat1/phat2
		phil = 0.95 * phihat
		r = 0
		if (x != 0) {
		  while (r <= zstar) {
			a = (n + m) * phil
			b = -((y + m) * phil + x + n)
			c = y + x
			p1hat = (-b - sqrt(b^2 - 4 * a * c))/(2 * a)
			p2hat = p1hat * phil
			q2hat = 1 - p2hat
			var = (n * m * p2hat)/(m * (phil - p2hat) + n * q2hat)
			r = ((x - m * p2hat)/q2hat)/sqrt(var)
			CIL = phil
			phil = CIL/1.0001
		  }
		}
		
		phiu = 1.1 * phihat
		if (x == 0) {
		  CIL = 0
		  if (n < 100) {
			phiu =  0.01
		  }
		  else {
			phiu = 0.001
		  }
		}
		
		r = 0
		while (r >= -zstar) {
		  a = (n + m) * phiu
		  b = -((y + m) * phiu + x + n)
		  c = y + x
		  p1hat = (-b - sqrt(b^2 - 4 * a * c))/(2 * a)
		  p2hat = p1hat * phiu
		  q2hat = 1 - p2hat
		  var = (n * m * p2hat)/(m * (phiu - p2hat) + n * q2hat)
		  r = ((x - m * p2hat)/q2hat)/sqrt(var)
		  phiu1 = phiu
		  phiu = 1.0001 * phiu1
		}
		
		CIU = phiu1
	  }
	}
	
	else if ((y != n) & (x != m) & (x != 0) & (y != 0)) {
	
		inits_l = rat*exp(-1*invnormal(1 - alpha/2)*sqrt(1/x + 1/y - 1/m - 1/n))
		inits_u = rat*exp(invnormal(1 - alpha/2)*sqrt(1/x  + 1/y - 1/n - 1/m))
		
		S = solvenl_init()
		solvenl_init_evaluator(S, &koopmancifun())
		solvenl_init_type(S, "zero")
		solvenl_init_technique(S, "newton")
		solvenl_init_numeq(S, 1)
		solvenl_init_narguments(S, 1)	
		solvenl_init_argument(S, 1, (x, m, y, n, alpha))
		solvenl_init_startingvals(S, inits_l)
		solvenl_init_iter_log(S, "off")
		CIL = solvenl_solve(S)
	  
		S = solvenl_init()
		solvenl_init_evaluator(S, &koopmancifun())
		solvenl_init_type(S, "zero")
		solvenl_init_technique(S, "newton")
		solvenl_init_numeq(S, 1)
		solvenl_init_narguments(S, 1)	
		solvenl_init_argument(S, 1, (x, m, y, n, alpha))
		solvenl_init_startingvals(S, inits_u)
		solvenl_init_iter_log(S, "off")
		CIU = solvenl_solve(S)
	}
}
	ci = (CIL, CIU)
	
	st_matrix("ci", ci)	
}
//************************************************************************	  
void U(real scalar a, real scalar r, real rowvector v){

	real scalar phat
	x = v[1] 
	m = v[2] 
	y = v[3] 
	n = v[4]
	alpha = v[5]
	
	phat = (a * (m + y) + x + n - ((a * (m + y) + x + n)^2 - 4 * a * (m + n) * (x + y))^0.5)/(2 * (m + n))
	r = (((x - m * phat)^2)/(m * phat * (1 - phat))) * (1 + (m * (a - phat))/(n * (1 - phat))) - invchi2(1, 1 - alpha)
}
//************************************************************************	
void U0(real scalar a, real scalar r, real rowvector v){

	m = v[2] 
	n = v[4]
	alpha = v[5]
	
	if (a <= 1) {
		r = m * (1 - a)/a - invchi2(1, 1 - alpha)
	}
	else {
		r = (n * (a - 1)) - invchi2(1, 1 - alpha)
	}
}
//************************************************************************	
/*
void numint_fe(real scalar len) {
	nu = st_data(., ("nu_b"))
	Ocoef = st_matrix("e(b)")
	ncols = cols(Ocoef)
	sigma_u = exp(Ocoef[ncols])

	pi = J(len, 1, .)
	u = range(0.001, 0.999, 0.001)
	z = invnormal(u)*sigma_u

	for (i=1; i<=len; i++) {
		nu_i = J(999, 1, nu[i])
		pi[i] = mean(invlogit(nu_i + z))
	}

	st_store(., "pi_int", pi)
}
*/
//************************************************************************
/*void numint_re(real scalar len) {
	nu = st_data(., ("nu_b"))
	Ocoef = st_matrix("e(b)")
	ncols = cols(Ocoef)
	sigma_u = exp(Ocoef[ncols-1])
	sigma_t = exp(Ocoef[ncols])

	pi = J(len, 1, .)
	u = range(0.001, 0.999, 0.001)
	z_u = invnormal(u)*sigma_u
	z_t = invnormal(u)*sigma_t

	for (i=1; i<=len; i++) {
		nu_i = J(999, 1, nu[i])
		pi[i] = mean(invlogit(nu_i + z_u + z_t))
	}

	st_store(., "pi_int", pi)
}*/

	mata mlib add lmetapreg *()
	mata mlib index
end
