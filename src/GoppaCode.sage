class GoppaCode:

	def __init__(self, n, m, g):
		t = g.degree();
		F2 = GF(2);
		F_2m = g.base_ring(); Z = F_2m.gen();
		PR_F_2m = g.parent(); X = PR_F_2m.gen();

		factor_list = list(factor(2^m-1)); final_factor_list = [];
		for i in range(len(factor_list)):
			for j in range(factor_list[i][1]):
				final_factor_list.append(factor_list[i][0]);

		while 1:
			primitive_root = F_2m.random_element();
			if primitive_root == 0:
				continue;
	
			for i in range(len(final_factor_list)):
				for j in itertools.combinations(final_factor_list,i):
					exponent = 1;			
					for _ in range(len(j)):
						exponent *= j[_];			
					if primitive_root^(exponent) == 1:
						output = False;
						break;
					else:
						output = True;
						continue;
				if output == False:
					break;
			if output == True:
				#print 'primitive root: ',primitive_root;
				break;
		#primitive_root = F_2m([0, 0, 1, 0, 1]);#TODO
		print "primitive_root", BinRepr(primitive_root);
		#Initialize the code locators
		codelocators = [];
		for i in range(2^m-1):
			codelocators.append(primitive_root^(i+1));
		codelocators.append(F_2m(0));
		print 'codelocators ', codelocators;
		#This is the same h to which Bernstein
		#refers in his polynomial view of a Goppa code
		h = PR_F_2m(1);
		#for a_i in codelocators:
			#h = h*(X-a_i);
		
		#Gamma is a list of the polynomials
		#used to determine membership in the code
		gamma = [];
		for a_i in codelocators:
			gamma.append((h*((X-a_i).inverse_mod(g))).mod(g));

		#Calculate the parity-check matrix (as polynomials)
		#In other words, the first column of this matrix
		#has entries that correspond to the coefficients of
		#the first polynomial of gamma.
		H_check_poly = matrix(F_2m, t, n);
		for i in range(n):
			#print '__init__: gamma[%d]=' %(i);
			#print gamma[i];
			#print '__init__: coeffs=';
			coeffs = list(gamma[i]);

			for j in range(t):
			# Check to make sure the coefficient exists
			#(It may not if the polynomial is not of
			#degree t)
				if j < len(coeffs):
					H_check_poly[j,i] = coeffs[j];
				else:
					H_check_poly[j,i] = F_2m(0);
		
		#Construct the binary parity-check matrix for the Goppa code.
		#Do so by converting each element of F_2^m to its binary
		#representation.
		H_Goppa = matrix(F2,m*H_check_poly.nrows(),H_check_poly.ncols());
		
		for i in range(H_check_poly.nrows()):
			for j in range(H_check_poly.ncols()):
				#print 'cnt=',i*n+j;
				#print '__init__: poly=',H_check_poly[i,j];
				be = bin(H_check_poly[i,j].integer_representation())[2:];
				be = be[::-1];
				#print '__init__: be=',H_check_poly[i,j].integer_representation();
				#be = '0'*(m-len(be))+be; 
				be = be+'0'*(m-len(be));
				be = list(be);
				H_Goppa[m*i:m*(i+1),j] = vector(map(int,be));
				#print "H_Goppa ",H_Goppa[m*i:m*(i+1),j].transpose();
		
                #Construct the generator matrix for our code by computing
		#a basis for the null-space of H_Goppa. The null-space is,
		#by definition, the codewords of our code.
		#G_Goppa = H_Goppa.right_kernel().basis_matrix();
		G_Goppa = H_Goppa.transpose().kernel().basis_matrix();
		G_Goppa_poly = H_check_poly.transpose().kernel().basis_matrix();                

		print 'k=n-mt= ',n-m*g.degree();
		print 'H_goppa_poly rank ', H_check_poly.rank();
		print H_check_poly.str();
		print 'H_goppa rank ', H_Goppa.rank();
		print H_Goppa.str();
		
		print 'G_goppa_poly nrows', G_Goppa_poly.nrows();	
		print G_Goppa_poly.str();
		print 'G_goppa nrows', G_Goppa.nrows();	
		print G_Goppa.str();
		#Construct the syndrome calculator. This will be used
		#to simplify the calculation of syndromes for decoding.
		SyndromeCalculator = matrix(PR_F_2m, 1, len(codelocators));
		for i in range(len(codelocators)):
			SyndromeCalculator[0,i] = (X - codelocators[i]).inverse_mod(g);

		#Remember these values
		self._n = n;
		self._m = m;
		self._g = g;
		self._t = t;
		self._codelocators = codelocators;
		self._SyndromeCalculator = SyndromeCalculator;
		self._H_Goppa = H_Goppa;
		self._H_gRS = H_check_poly;
		self._G_Goppa = G_Goppa;

		self._R = 0;self._alpha = 0;self._beta = 0;self._sigma = 0;
		
	def Encode(self, message):
		#Encoding a k-bit binary message done by
		#multiplication on the right by the generator matrix.
		return (message*self._G_Goppa);

	def _split(self,p):
		# split polynomial p over F into even part po
		# and odd part p1 such that p(z) = p2 (z) + z p2 (z)
		Phi = p.parent()
		p0 = Phi([sqrt(c) for c in p.list()[0::2]]);
		p1 = Phi([sqrt(c) for c in p.list()[1::2]]);
		return (p0,p1);

	def _g_inverse(self, p):
		# returns the g-inverse of polynomial p
		(d,u,v) = xgcd(p,self.goppa_polynomial());
		return u.mod(self.goppa_polynomial());

	def _norm(self,a,b):
	#This is the way in which Bernstein indicates
	#the norm of a member of the lattice is
	#to be defined.
		X = self.goppa_polynomial().parent().gen();
		return 2^((a^2+X*b^2).degree());

	def _lattice_basis_reduce(self, s):
		# a <- s   b <- v 
		g = self.goppa_polynomial();
		t = g.degree();
		a = []; a.append(0);
		b = []; b.append(0);
		(q,r) = g.quo_rem(s);
		(a[0],b[0]) = simplify((g - q*s, 0 - q))
			
			
		#print 'lattice: a[0]=',BinRepr(a[0]);
		#print 'lattice: b[0]=',BinRepr(b[0]);
		

		#If the norm is already small enough, we
		#are done. Otherwise, intialize the base
		#case of the recursive process.
		#print 'norm: ',self._norm(a[0],b[0]);
		if self._norm(a[0],b[0]) > 2^t:
			a.append(0); b.append(0);
			(q,r) = s.quo_rem(a[0]);
			(a[1],b[1]) = (r, 1 - q*b[0]);
			if a[1] == 0:
				return (s,1);	
			#print 'lattice: a[1]=',BinRepr(a[1]);
			#print 'lattice: b[1]=',BinRepr(b[1]);
			
		else:
			return (a[0], b[0]);
		#Continue subtracting integer multiples of
		#the shorter vector from the longer until
		#the produced vector has a small enough norm.
			
		i = 1;
		while self._norm(a[i],b[i]) > 2^t:
			a.append(0); b.append(0);
			(q,r) = a[i-1].quo_rem(a[i]);
			(a[i+1],b[i+1]) = (r, b[i-1] - q*b[i]);
			i+=1;
			#print 'lattice: a[',i,']= ',BinRepr(a[i]);
			#print 'lattice: b[',i,']= ',BinRepr(b[i]);
		return (a[i],b[i]);


	
	def _extended_euclidean(self, a,b,degree):
		# v*b = s mod a, a>b
		s = []; u = []; v = [];
		s.append(a); s.append(b);
		u.append(1); u.append(0);
		v.append(0); v.append(1);
		i=1;
		
		while s[i].degree() >= degree:
			i += 1;	
			s.append(0);u.append(0);v.append(0);
			(q,r) = s[i-2].quo_rem(s[i-1]);		
			s[i] = s[i-2] + q*s[i-1];
			u[i] = u[i-2] + q*u[i-1];
			v[i] = v[i-2] + q*v[i-1];
			#print 'euclidean: omega[i]=',BinRepr(s[i]);
			#print 'euclidean: sigma[i]=',BinRepr(v[i]);
		sigma = v[i]
		omega = s[i]
		return (sigma, omega);	 
	
	def _extended_euclidean_I(self,a,b,degree):
		s = []; u = []; v = [];
		s.append(a); s.append(b);
		u.append(1); u.append(0);
		v.append(0); v.append(1);
		i=1;	

		while s[i].degree() >= degree:
			i += 1;
			q = 0;
			t = s[i-2].degree()-s[i-1].degree();
			s.append(0);u.append(0);v.append(0);
			while t >= 0:
				coeff1 = list(s[i-2])[-1];
				coeff2 = list(s[i-1])[-1];
				X = a.parent().gen();
				Q = coeff1 * 1/coeff2 * X^t;
				s[i] = s[i-2] + Q*s[i-1];
				q += Q;  
				s[i-2] = s[i];
				t = s[i-2].degree()-s[i-1].degree(); 
			u[i] = u[i-2] + q*u[i-1];
			v[i] = v[i-2] + q*v[i-1];

		sigma = v[i]
		omega = s[i]
		return (sigma, omega);
	
	def _extended_euclidean_II(self,a,b,degree):
		s = []; u = []; v = [];
		s.append(a); s.append(b);
		u.append(1); u.append(0);
		v.append(0); v.append(1);
		i=1;	
		
		while s[i].degree() >= degree:
			i += 1;
			delta = 1; q = 0;
			t = s[i-2].degree()-s[i-1].degree();
			s.append(0);u.append(0);v.append(0);
			tmp = s[i-2];
			while t >= 0:
				coeff1 = list(s[i-2])[-1];
				coeff2 = list(s[i-1])[-1];
				X = a.parent().gen();
				Q = coeff1 * X^t;
				s[i] = coeff2*s[i-2] - Q*s[i-1];
				delta *= coeff2;				
				q = coeff2*q + Q;  
				s[i-2] = s[i];
				t = s[i-2].degree()-s[i-1].degree(); 
			assert s[i] == delta*tmp - q*s[i-1], 'relation wrong!';
			u[i] = delta*u[i-2] - q*u[i-1];
			v[i] = delta*v[i-2] - q*v[i-1];

		sigma = v[i];
		omega = s[i];
		return (sigma, omega);	

		
	def Decode(self, word_, mode='Patterson'):
		g = self._g;
		word = copy(word_);
		X = g.parent().gen();

		#Compute the syndrome necessary for Patterson’s Algorithm.
		synd = self._SyndromeCalculator*word.transpose();		
		#print 'Decode: synd=', synd;
	
		syndrome_poly = 0;
		for i in range (synd.nrows()):
			syndrome_poly += synd[i,0]*X^i
		
		if mode == 'Patterson':
			#We will decode codewords using Patterson’s Algorithm, error-locator polynomial returned.
			error = self.SyndromeDecode(syndrome_poly, 'Patterson');
			return word+error;
		
		elif mode == 'Euclidean':
			error = self.SyndromeDecode(syndrome_poly, 'Euclidean');	
			return word+error;
			
	
	def SyndromeDecode(self, syndrome_poly, mode='Patterson'):
		#Unlike Decode(), SyndromeDecode uses the syndrome polynomial as its input, outputs the error pattern.
			
		g = self.goppa_polynomial();
		X = g.parent().gen();
		error = matrix(GF(2),1,self.parity_check_matrix().ncols());
		
		if mode == 'Patterson':	
			#Compute the syndrome necessary for Patterson’s Algorithm.
			#Take the necessary square root
			(g0,g1) = self._split(g); sqrt_X = g0*self._g_inverse(g1);
		        T = syndrome_poly.inverse_mod(g);
		        		
			(T0,T1) = self._split(T - X);
			R = (T0+ sqrt_X*T1).mod(g);
		
			#Perform lattice basis reduction.
			(alpha, beta) = self._lattice_basis_reduce(R);

			#Construct the error-locator polynomial.
			sigma = (alpha*alpha) + (beta*beta)*X;
			print 'SyndromeDecode: sigma=', BinRepr(sigma);
			#Pre-test sigma
			if (X^(2^m)).mod(sigma) != X:	
				print "sigma: Decodability Test Failed";
				return error;# return a zero vector 
			#For every root of the error polynomial,
			#correct the error induced at the corresponding index.
			for i in range(len(self._codelocators)):
				if sigma(self._codelocators[i]) == 0:
					error[0,i] = 1;
			return error;

		elif mode == 'Euclidean':
			(sigma, omega) = self._extended_euclidean_II(g,syndrome_poly,floor(g.degree()/2));
			print "Decode, sigma=",sigma;
			sigma = sigma/sigma.coefficients()[-1];
			#For every root of the error polynomial,
			#correct the error induced at the corresponding index.
			for i in range(len(self._codelocators)):
				if sigma(self._codelocators[i]) == 0:
					error[0,i] = 1;
			return error;
			
	#Accessors
	def generator_matrix(self):
		return (self._G_Goppa);

	def goppa_polynomial(self):
		return (self._g);

	def parity_check_matrix(self):
		return (self._H_Goppa);
	
	def parity_check_poly_matrix(self):
		return (self._H_gRS);
	
	def R(self):
		return self._R;
	def error_locator(self):
		return (self._sigma, self._alpha, self._beta);

	


