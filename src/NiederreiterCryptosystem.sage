class NiederreiterCryptosystem:

	def __init__(self, n, m, g, delta):
		#Construct the Goppa code
		goppa_code = GoppaCode(n,m,g);
		
		#assert goppa_code.generator_matrix().nrows() == n-m*g.degree(), "Generator Matrix total row number is incorrect.";
		k = goppa_code.generator_matrix().nrows();
		print 'genertor.nrows',k;
		# Set up the random scrambler matrix
		S = matrix(GF(2),n-k,[random()<0.5 for _ in range((n-k)^2)]);
		while (rank(S) < n-k) :
			S[floor((n-k)*random()),floor((n-k)*random())] +=1;
		
		# Set up the permutation matrix
		rng = range(n); P = matrix(GF(2),n);
		for i in range(n):
			p = floor(len(rng)*random());
			P[i,rng[p]] = 1; rng=rng[:p]+rng[p+1:];
		
		#Remember these values
		self._m_GoppaCode = goppa_code;
		self._g = g;
		self._t = g.degree();
		self._S = S;
		self._P = P;
		self._PublicKey = S*(self._m_GoppaCode.parity_check_matrix())*P;
		self._delta = delta;
		
	#This is a help function which will be useful for encryption.
	def _GetRowVectorWeight(self,n):
		weight = 0;
		for i in range(n.ncols()):
			if n[0,i] == 1:
				weight = weight+1;
		return weight;

	
	def Encrypt(self, message):
	#Assert that the message is of the right length
		assert (message.ncols() == self.public_key().ncols()), "Message is not of the correct length."; 
		#assert (self._GetRowVectorWeight(message) <= self.max_num_errors()), "Message weighs at most %d." %(self.max_num_errors());
		
		code_word = self.public_key()*(message.transpose());
		
		return code_word.transpose();
	
	def Decrypt(self, received_word):
	#Assert that the received_word is of the right length
		received_word = received_word.transpose();
		assert (received_word.nrows() == self.public_key().nrows()), "Received word is not of the correct row length.";
                assert (received_word.ncols() == 1), "Received word is not of the correct column length."
		#Strip off the random scrambler matrix and decode the received word.
		message = ~(self._S)*received_word;
		
		
		#Convert message into syndrome polynomial.
		t = self.max_num_errors();
		m = message.nrows()/t;
		g = self.goppa_code().goppa_polynomial();		
		F2 = GF(2);
		F_2m = g.base_ring(); Z = F_2m.gen();
		PR_F_2m = g.parent(); X = PR_F_2m.gen();		
		syndrome_poly = 0;
                
		for i in range(t):
                        tmp = [];
			for j in range(m):
				tmp.append(message[i*m+j,0]);
                        syndrome_poly += F_2m(tmp[::1])*X^i;
               
		#Retrieve using syndrome decoding algorithm. 		
		message = self.goppa_code().SyndromeDecode(syndrome_poly);
		
		#Solve the system to determine the original message.
		message = message*self._P;
		return message;			

	
	def _str2num(self,ascii_text):
		#Turn ascii text into big integers by treating the ASCII values 
		#of each character as digits in base 256:
    		return ZZ(map(ord,ascii_text),256);

	def _common_primroot(self,p,q):
	    # Finds a primitive root common to both p and q
	    ps = [f for f,e in factor(p-1)]
	    qs = [f for f,e in factor(q-1)]
	    a=2
	    while true:
		outp = true
		for d in ps:
		    if power_mod(a,(p-1)//d,p)==1:
		        outp=false
		outq =  true
		for d in qs:
		    if power_mod(a,(q-1)//d,q)==1:
		        outp=false
		if outp and outq:
		    return a
		    break
		else:
		    a+=1 

	def CryptoHash(self, plaintext, hash_value_length):
		#This cryptographic hash function was proposed by Adi Shamir.
		#For more details, visit http://amca01.wordpress.com/2010/09/03/a-quick-and-dirty-hash-function/
		# Implements Shamir's hash function h = a^m (mod p*q) where a is a
    		# number of high multiplicative order mod p*q.  This works well if
    		# a is chosen to be a common primitive root of both p and q
		p = next_prime(2^150);
		q = next_prime(2^151);
		a = self._common_primroot(p,q);
		
		hashvalue = power_mod(a,self._str2num(plaintext),p*q)%(2^hash_value_length);
		hashvalue_txt = bin(hashvalue)[2:];
		hashvalue_txt = '0'*(hash_value_length-len(hashvalue_txt))+hashvalue_txt; 
    		return matrix(map(GF(2),hashvalue_txt));
		
	def Signature(self, message, mode='random'):
		#decrypt a message of length m*t using complete docoding algorithm, here we try to decode a 't+delta' code.
		new_message = matrix(GF(2), 1, message.ncols());
		cnt = 0;
		if mode == 'random':
			#randomly select columns from PublicKey Matrix			
			while 1:
				print 'Signature:',cnt;
				cnt += 1;
				new_message = message;	random_pos_list = [];
				rng = range(message.ncols());
					
				for i in range(self._delta):
					p = floor(len(rng)*random());
					new_message += self._PublicKey[:,rng[p]].transpose(); 
					rng=rng[:p]+rng[p+1:]; random_pos_list.append(p);
			
				new_message = self.Decrypt(new_message);
			
				if self._GetRowVectorWeight(new_message) == self._t:
					# a decodable hash value is found.
                                        for i in random_pos_list:
						new_message[0,i] += 1;
                                        if self._GetRowVectorWeight(new_message) == self._t+self._delta:
                                                print 'found:', new_message.str(), random_pos_list;
                                                break;
                                        else:
                                                print 'not found, because weight is not desired.';
				else: 
					print 'not found, weight:',self._GetRowVectorWeight(new_message);
		
			return new_message;
		
		elif mode == 'sequential':
			#Sequentially select columns from PublicKey Matrix
			rng = range(self.public_key().ncols());
			for i in itertools.combinations(rng,self._delta):
				print 'Signature:',cnt;	cnt += 1;
				new_message = message;	random_pos_list = [];
				
				for j in range(self._delta):				
					new_message += self._PublicKey[:,i[j]].transpose();
					random_pos_list.append(i[j]);
				
				new_message = self.Decrypt(new_message);		
				if self._GetRowVectorWeight(new_message) == self._t:
				# a decodable hash value is found.
					for i in random_pos_list:
						new_message[0,i] += 1; 
					print 'found:', new_message.str(), random_pos_list;
					return new_message;
				else: 
					print 'not found, weight:', self._GetRowVectorWeight(new_message);
						
			#could not decode this hashvalue					
			print 'Signature: failed!'; return None;
		
		else:
			print 'Signature(self, message, mode=''random''): Invalid mode'; pass;
		
		
	
	def Verify(self, signature, ascii_text_hash):
		#Encrypt signature and then compare it with the document's hash value.
		if self.Encrypt(signature) == ascii_text_hash :
			#This is a valid signature.
			return True;
		else :
			#This is an invalid signature.			
			return False;			

	#Accessors
	def public_key(self):
			return self._PublicKey;

	def goppa_code(self):
			return self._m_GoppaCode;

	def max_num_errors(self):
			return self._t;
		
