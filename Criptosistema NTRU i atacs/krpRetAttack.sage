load("auxMethods.sage")

import logging
import random
import numpy as np
import itertools
import time
import base64
import sympy

class Krp_lattice:

    N = None
    p = None
    q = None
    d = None
    r = None
    h = None
    f = None
    f_inv_mod_p = None
    f_inv_mod_q = None
    g = None
    r_ring = None
    r_ring_mod_p = None
    r_ring_mod_q = None
    r_ring_mod_2 = None

    def __init__(self, filename): #Inicialitzacio dels parametres i dels anells que es requereixen per fer funcionar
                                  #l'atac al KRP basat en reticles.

        f = open(filename, "r") #Es llegeix el fitxer que conte els parametres publics.
        params = f.read()
        f.close()

        param_list = [] #Llista que conte els parametres publics en format string.

        param_list  = params.split("A")

        param_list[-1] = param_list[-1].split(",") #Es separen els coeficients de la clau publica h llegida del fitxer.

        h_int_list = [int(c) for c in param_list[-1]] #Es passen de coeficients de tipus string a tipus int.
        z_pol_ring.<x> = PolynomialRing(ZZ)

        self.N = Integer(param_list[0]) #S'estableixen els nous valors dels parametres publics.
        self.p = Integer(param_list[1])
        self.q = Integer(param_list[2])
        self.d = Integer(param_list[3])
        self.h = z_pol_ring(h_int_list)

        aux_ring.<x> = PolynomialRing(ZZ)
        self.r_ring = aux_ring.quotient(x^self.N - 1)

        aux_ring.<x> = PolynomialRing(ZZ.quotient(self.p)) #PolynomialRing(QuotientRing(ZZ, ZZ.ideal(self.p)))
        self.r_ring_mod_p = aux_ring.quotient(x^self.N - 1) #QuotientRing(aux_ring, aux_ring.ideal(x^self.N - 1))

        aux_ring.<x> = PolynomialRing(ZZ.quotient(self.q)) #PolynomialRing(QuotientRing(ZZ, ZZ.ideal(self.q)))
        self.r_ring_mod_q = aux_ring.quotient(x^self.N - 1) #QuotientRing(aux_ring, aux_ring.ideal(x^self.N - 1))

        aux_ring.<x> = PolynomialRing(ZZ.quotient(2))
        aux_ring.<x> = PolynomialRing(ZZ.quotient(2))
        self.r_ring_mod_2 = aux_ring.quotient(x^self.N - 1)



    def decryption(self, e): #Es desencripta el missatge entrat per parametre.

        z_pol_ring.<x> = PolynomialRing(ZZ)
        a = z_pol_ring(list(self.f))*e             #Es calcula la primera de les dues convolucions.
        a = center_lift(self.r_ring_mod_q(a)) #Es prenen moduls q i (x^N-1) i es centre liften els coeficients.
        b = z_pol_ring(list(a))*self.f_inv_mod_p   #Es calcula la segona de les dues convolucions.
        b = center_lift(self.r_ring_mod_p(b)) #Es prenen moduls p i (x^N-1) i es centre liften els coeficients.

        return b



    def decrypt_message(self, list_of_encrypted): #Es desencripten el conjunt de submissatges encriptats.

        decrypted_list = []

        for i in range(len(list(list_of_encrypted))): #Per cadascun dels submissatges es crida a la funcio que els
                                                      #desencripta.

            decrypted_list.append(list(self.decryption(list_of_encrypted[i])))

        return decrypted_list #Es retorna una llista amb tots els submissatges desencriptats.


    def decrypt_file(self, e_filename, d_filename): #Donats dos fitxers, llegeix els coeficients dels polinomis
                                                    #corresponents al text encriptat i escriu el missatge desencriptat al
                                                    #segon fitxer.

        decrypted_lines_list = [] #Llista on els elements son les linies de text desencriptades.
        b_list = [] #Llista que conte subllistes corresponents als bits dels submissatges d'una linia de text
                    #desencriptada.
        b_def_list = [] #Llista que conte els bits dels submissatges unificats d'una linia de text desencriptada.
        pol_lines_list = [] #Llista que conte subllistes referents als coeficients dels missatges encriptats per cada linia
                            #de text.
        pol_lines_list = read_encrypted_file(e_filename, self.N) #Es llegeix el text encriptat del fitxer corresponent.

        for line in pol_lines_list: #Es desencripten les linies del fitxer una a una.
            b_list = a.decrypt_message(line) #Es desencripta una linia.
            add_zeros(b_list, self.N) #S'afegeixen els 0 que calgui al final de cada polinomi.
            b_def_list = [bit for polynomial in b_list for bit in polynomial] #S'unifiquen els coeficients dels polinomis
                                                                              #a una sola linia.
            dec_text = "".join(map(chr, np.packbits(np.array(list(b_def_list))))) #Es transformen els bits al text en
                                                                                  #string corresponent.
            decrypted_lines_list.append(dec_text+"\n") #S'afegeix un salt de linia al final de cada linia.

        f = open(d_filename, "w") #S'escriu el text desencriptat al fitxer corresponent.
        f.writelines(decrypted_lines_list)
        f.close()



    def attack_file(self, e_filename, a_filename): #S'ataca el primer fitxer entrat per parametre i es guarda el missatge
                                                   #al segon fitxer entrat per parametre.

        self.obtain_priv_keys() #S'obtenen les claus privades mitjançant l'atac.

        self.decrypt_file(e_filename, a_filename) #Es desencripta el fitxer.


    def obtain_priv_keys(self): #S'obtenen les claus privades seguint l'atac a KRP basat en reticles.

        keys_pair_list = [] #Llista que conte els coeficients de les claus privades f i g.
        f_list = [] #Llista que conte els coeficients de la clau privada f.
        g_list = [] #Llista que conte els coeficients de la clau privada g.

        M = self.gen_M_matrix() #Es genera la matriu M que defineix el reticle.
        reduced_M = M.LLL(algrithm='NTL:LLL', fp= None) #S'aplica l'algorisme LLL per reduir la base del reticle.

        keys_pair_list = list(reduced_M[0]) #S'obtenen les claus privades f i g.


        f_list = keys_pair_list[:self.N] #Es separen les claus privades f i f.
        g_list = keys_pair_list[self.N:]

        #print(list(self.h))
        #print(M)
        #print("\n")
        #print("AAAAAAAAAAAAAAAA")
        #print("\n")
        #print(reduced_M)
        #print(keys_pair_list)
        #print(f_list)
        #print(g_list)

        self.f = self.r_ring(f_list) #S'emmagatemen les claus privades f i g i els inversos modul p i q de f.
        self.g = self.r_ring(g_list)

        self.f_inv_mod_q = inv_pol_mod_power_2(self.r_ring_mod_q(f_list), self.r_ring_mod_2(f_list), self.N, self.q)

        if not self.p.is_power_of(2):

            self.f_inv_mod_p = inv_pol_mod_prime(self.r_ring_mod_p(f_list))

        else:
            self.f_inv_mod_p = inv_pol_mod_power_2(self.r_ring_mod_p(f_list), self.r_ring_mod_2(f_list), self.N, self.p)


    def gen_M_matrix(self): #Es genera la matriu M que defineix el reticle.

        Id = matrix.identity(self.N) #Es genera la submatriu identitat.
        q_matrix = self.q*Id #Es genera la submatriu q*identitat.
        zero_matrix = matrix(self.N) #Es genera la submatriu de zeros.

        h_matrix = self.gen_h_matrix() #Es genera la submatriu de la clau publica h.

        M = block_matrix(ZZ,[[Id,h_matrix],[zero_matrix,q_matrix]]) #S'ajunten les 4 matrius per construir la matriu M.

        return M

    def gen_h_matrix(self): #Es genera la submatriu de la clau publica h.

        h_rows_list = [] #Llista que conte totes les files de la submatriu.

        h_list = list(self.h)

        if len(h_list) < self.N: #Si la llista dels seus coeficients es menor que N, s'afegeixen els 0 que calgui al
                                 #final.
            zeros = [0]*(self.N - len(h_list))
            h_list.extend(zeros)

        h_rows_list.append(h_list) #S'afegeix la primera fila de la matriu, que correspon al polinomi h

        for i in range(1,self.N): #S'itera per la resta de files

            h_list = h_list[-1:] +  h_list[:-1] #Es roten en una posicio cap a la dreta els coeficients del polinomi h.
            h_rows_list.append(h_list) #S'afegeix la seguent fila de la matriu.

        h_matrix = matrix(ZZ, h_rows_list) ##Es genera la submatriu de la clau publica h a partir de les files creades.

        return h_matrix





