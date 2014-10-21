#!/home/david/sage-6.2-i686-Linux/sage -python

import pdb
import itertools
import time
from sage.all import *

############################################
#Functions & Classes to Import
############################################
load("./GoppaCode.sage")
load("./McElieceCryptosystem.sage")
load("./NiederreiterCryptosystem.sage")

def GetRandomMessage(message_length):
	message = matrix(GF(2), 1, message_length);
	for i in range(message_length):
		message[0,i] = randint(0,1);
	return message;

def GetRandomMessageWithWeight(message_length, message_weight):
	message = matrix(GF(2), 1, message_length);
	rng = range(message_length);
	for i in range(message_weight):
		p = floor(len(rng)*random());
		message[0,rng[p]] = 1; rng=rng[:p]+rng[p+1:];
	return message;

def GetGoppaPolynomial(polynomial_ring, polynomial_degree):
	while 1:
		irr_poly = polynomial_ring.random_element(polynomial_degree);
		irr_poly_list = irr_poly.list(); irr_poly_list[-1] = 1;
		irr_poly = polynomial_ring(irr_poly_list);
		if irr_poly.degree() != polynomial_degree:
			continue;
		elif irr_poly.is_irreducible():
			break;				
		else :
			continue;
	
	return irr_poly;

def BinRepr(poly):
	try: 
		poly_ls = list(poly);
	except TypeError:
		bin_repr = bin(poly.integer_representation())[2:];
		bin_repr = bin_repr[::-1];
		print bin_repr;
		return;
		
	else: 
		for _ in poly_ls:
			bin_repr = bin(_.integer_representation())[2:];
			bin_repr = bin_repr[::-1];
			print bin_repr,
			
	
############################################
#Beginning of Script
############################################
#The goal of the script will be the testbech for: 
#1) GoppaCode Encoding&Decoding
#2) McElieceCryptosystem Encryption&Decryption
#3) NiedereiterCryptosystem Encryption&Decryption and Signature&Verification

#specifiy your testing option
option = 'GoppaCode'; 

#pdb.set_trace();
for k in range(1):
	
	# setup system parameters (m,n,t, goppa_polynomial)
	m = 4;
	n = 2**m;
	t = 2;
	delta = 2;
	#construct finite field and its extension[[1 1 0 1 1 0 0 0 1]]
	#P = PolynomialRing(GF(2),'M'); x = P.gen(); f = 1+x+x**3+x**4+x**8;#TODO
	F_2m = GF(n,'Z', modulus='random'); print 'modulus=', F_2m.modulus(); 
	PR_F_2m = PolynomialRing(F_2m,'X');
	Z = F_2m.gen(); X = PR_F_2m.gen(); 
		
	for _ in range(1):	
		irr_poly = GetGoppaPolynomial(PR_F_2m, t);
		#irr_poly = PR_F_2m([F_2m([1,0,1,1,0,0,1]),F_2m([1,0,0,0,1,0,1,1]),F_2m([0,0,1,0,0,0,1,1]),F_2m([1,1,0,0,0,1,0,1]),F_2m([1,0,1,0,0,1]),F_2m([0,1,1,1,1,0,1,1]),F_2m([1,1,0,1,0,1,1]),F_2m([1,1,1,1,1,0,0,1]),F_2m([0,1,0,0,0,1,0,1]),F_2m([0,1,0,0,0,1,1,1]),F_2m([0,1,0,1,1,1,1]),F_2m([1,0,0,0,1,1,1,1]),F_2m([0,1,1]),F_2m([0,0,1,0,1,0,1,1]),F_2m([1,0,0,1,1,0,1]),F_2m([0,0,0,1,1,0,1,1]),F_2m([0,1,1,1,1,1]),F_2m([0,0,0,1,0,0,1]),F_2m([0,0,1,1,0,0,1]),F_2m([0,1,0,0,1]),F_2m([1,1,0,0,1]),F_2m([1,1,1,1,1,0,1]),F_2m([0,0,0,1]),F_2m([1,0,1,1,0,0,0,1]),F_2m([1,0,1,0,0,0,0,1]),F_2m([0,1,1,0,1]),F_2m([1,0,0,0,0,0,0,1]),F_2m([1,0,0,1,1,0,1]),F_2m([1,1,0,1,1,1]),F_2m([1,0,1,1,0,1,1,1]),F_2m([1,0,1,1,0,1,1,1]),F_2m([1])]);#TODO
		#irr_poly = irr_poly*irr_poly;
		#t = t*2;		
		print 'm=%d, n=%d, t=%d' %(m,n,t);		
		print 'Goppa-polynomial:',irr_poly;

		if option == 'GoppaCode':
			goppacode = GoppaCode(n,m,irr_poly);
			#Decoding GoppaCode using classic matrix form.
			#Get a random message and encrypt it.
			message = GetRandomMessage(goppacode.generator_matrix().nrows());
			#message = matrix(GF(2),[0, 1, 0, 0, 1, 0, 0, 0]);#TODO
			codeword = goppacode.Encode(message);
			error = GetRandomMessageWithWeight(goppacode.generator_matrix().ncols(),t/2);
			#TODO                        
			#error = matrix(GF(2),[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1]);
			ciphertext =  codeword + error;

			print 'random message is:', message; print 'codeword is:', codeword.str();
			print 'error is:', error.str();
			print 'ciphertext is:', ciphertext.str();
			
			#Decrypt the ciphertext using algebraic syndrome decoding algorihtm 
			recoveredtext1 = goppacode.Decode(ciphertext,'Euclidean');
			
			#Decrypt the ciphertext using classic syndrome decoding algorithm
			g = goppacode.goppa_polynomial();	X = g.parent().gen();		
			syndrome = goppacode.parity_check_matrix()*ciphertext.transpose();		
			syndrome_poly = 0;		
			for i in range(t):
		                tmp = [];
				for j in range(m):
					tmp.append(syndrome[i*m+j,0]);
		                
				syndrome_poly += F_2m(tmp[::1])*(X**i);
			print 'syndrome_poly=', BinRepr(syndrome_poly);	
			recoveredtext2 = goppacode.SyndromeDecode(syndrome_poly,'Euclidean');
			
			sigma = 1;			
			for i in range(error.ncols()):
				if (error[0,i] == 1):
					sigma = sigma*(X-goppacode._codelocators[i]);
			print 'sigma= ',sigma;
			omega = 0; 
			for i in range(error.ncols()):
				if (error[0,i] == 1):
					tmp = sigma/(X-goppacode._codelocators[i]);
					omega = omega+tmp;
					

			print 'recovered message1 is:', recoveredtext1.str();
			print 'recovered message2 is:', recoveredtext2.str();
			
			if recoveredtext1 == codeword and recoveredtext2+ciphertext == codeword:
				print 'It works!';
			else :
				print 'Something wrong!'; 
				for i in range(len(goppacode._codelocators)):
					if goppacode._codelocators[i] == 1:
						print 'codelocators',i+1;
				raw_input('look at here!');
	
		elif option == 'McElieceCryptosystem':
			crypto = McElieceCryptosystem(n,m,irr_poly);
	
			#Get a random message and encrypt it.
			message = GetRandomMessage(crypto.goppa_code().generator_matrix().nrows());
			encrypted_message = crypto.Encrypt(message);

			#Decrypt the secret message
			decrypted_message = crypto.Decrypt(encrypted_message);

			print 'random msg:', message;
			print 'encrypted msg:', encrypted_message;
			print 'decrypted_msg:', decrypted_message;
			if message == decrypted_message:
				print 'It works!';
			else :
				print 'Something wrong!';
				raw_input('look at here!');

		
		elif option == 'NiedereiterCryptosystem':		
			crypto = NiederreiterCryptosystem(n,m,irr_poly,delta);
			#Encrypt & Decrypt			
			#Get m*t-bits random message weighing at most t and encrypt it.
			message = GetRandomMessageWithWeight(crypto.public_key().ncols(),irr_poly.degree());
			encrypted_message = crypto.Encrypt(message);
		        decrypted_message = crypto.Decrypt(encrypted_message);
			print 'random msg:', message.str();
			print 'encrypted msg:', encrypted_message.str();
		        print 'decrpted msg:', decrypted_message.str();
			if message == decrypted_message:
				print 'It works!';
			else :
				print 'Something wrong!', crypto._GetRowVectorWeight(decrypted_message);

			#Signature & Verification
			ascii_text = 'Hello, world! Hello, David!';
			hashval = crypto.CryptoHash(ascii_text, m*t);
			print 'text=', ascii_text;
			print 'hash=', hashval;
			
			sig = crypto.Signature(hashval,'random');
			if sig == None:
				# signature failure
                                print 'Signature failure!';
				exit();
			print 'sig=', crypto._GetRowVectorWeight(sig);
			if crypto.Verify(sig,hashval) == True:
				print 'Valid!';
			else:
				print 'Invalid!';

		else : 
			print 'The input option is incorrect!';

		
		
 
		
	

