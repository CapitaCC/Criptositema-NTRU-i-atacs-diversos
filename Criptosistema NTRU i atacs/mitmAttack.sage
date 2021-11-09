load("auxMethods.sage")

import logging
import random
import numpy as np
import itertools
import time
import base64
import sympy


class Meet_in_the_middle:

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
                                  #l'atac per trobada a mig cami.

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


    def obtain_priv_keys(self): #S'obtenen les claus privades seguint l'atac per trobada a mig cami.

        g_list =[0]*self.N #Llista que conte els coeficients de la clau privada g.

        a1_N = ceil(self.N/2) #Llista que conte els coeficients de la primera part de la clau privada f.
        a2_N = floor(self.N/2) #Llista que conte els coeficients de la segona part de la clau privada f.

        if self.d%2 == 0: #Repartiment dels 1 i -1 entre les parts del polinomi f quan d es parell.
            a1_one_list = [1]*Integer(self.d/2+1)
            a1_neg_one_list = [-1]*Integer(self.d/2)
            a1_zero_list = [0]*Integer(a1_N-self.d-1)
            a1_list = a1_one_list + a1_neg_one_list + a1_zero_list #Es concatenen el conjunt de 1, -1 i 0.
            a2_one_list = [1]*Integer(self.d/2)
            a2_neg_one_list = [-1]*Integer(self.d/2)
            a2_zero_list = [0]*Integer(a2_N-self.d)
            a2_list = a2_one_list + a2_neg_one_list + a2_zero_list #Es concatenen el conjunt de 1, -1 i 0.

        else: #Repartiment dels 1 i -1 entre les parts del polinomi f quan d es senar.
            a1_one_list = [1]*Integer(Integer(self.d+1)/2)
            a1_neg_one_list = [-1]*Integer(Integer(self.d+1)/2)
            a1_zero_list = [0]*Integer(a1_N-self.d-1)
            a1_list = a1_one_list + a1_neg_one_list + a1_zero_list #Es concatenen el conjunt de 1, -1 i 0.
            a2_one_list = [1]*Integer(Integer(self.d+1)/2)
            a2_neg_one_list = [-1]*Integer(Integer(self.d-1)/2)
            a2_zero_list = [0]*Integer(a2_N-self.d)
            a2_list = a2_one_list + a2_neg_one_list + a2_zero_list #Es concatenen el conjunt de 1, -1 i 0.

        z_pol_ring.<x> = PolynomialRing(ZZ)

        #a1_p = sorted(sage.combinat.permutation.Permutations_mset(a1_list)) #S'obtenen totes les possibles permutacions
                                                                            #de la primera part de la clau privada f.
        #a2_p = sorted(sage.combinat.permutation.Permutations_mset(a2_list)) #S'obtenen totes les possibles permutacions
                                                                            #de la segona part de la clau privada f.
        #a1_iter = itertools.permutations(a1_list)
        #a2_iter = itertools.permutations(a2_list)

        a1_iter = sympy.utilities.iterables.multiset_permutations(a1_list)
        a2_iter = sympy.utilities.iterables.multiset_permutations(a2_list)

        k = ceil(log(binomial(self.N/2,self.d/2),2)) #Es calcula el nombre de calaixos que s'usaran per classificar.
        bin_pol = [] #Llista que conte tots els polinomis d'un mateix calaix.
        bins = [bin_pol]*pow(2,k) #Llista de tots els calaixos.

        #while a1_p: #Es classifiquen totes les possibles primeres parts del polinomi f.

        for a1_p in a1_iter:

            #pol_a1_list = list(a1_p.pop())

            pol_a1_list = list(a1_p)
            pol_a1 = z_pol_ring(pol_a1_list)
            pol_g1 = pol_a1*self.h  #Es calcula la convolucio.
            pol_g1 = self.r_ring_mod_q(pol_g1).lift().change_ring(z_pol_ring) #Es prenen moduls p i (x^N-1)
                                                                              #i es centre liften els coeficients.
            g1_list = list(pol_g1)

            if len(g1_list) < self.N: #Si la llista dels seus coeficients es menor que N, s'afegeixen els 0 que calgui al
                                      #final.
                zeros = [0]*(self.N - len(g1_list))
                g1_list.extend(zeros)

            g1_k_coef_list = g1_list[:k] #Es prenen els primers k coeficients del polinomi g1.
            self.classify(pol_a1_list, g1_k_coef_list, bins) #Es classifica el polinomi g1.

        #while a2_p: #Es classifiquen totes les possibles segones parts del polinomi f

        for a2_p in a2_iter:

            #pol_a2_list = list(a2_p.pop())

            pol_a2_list = list(a2_p)
            neg_pol_a2_list = [ -x for x in pol_a2_list] #Es converteix el polinomi en negatiu.
            pol_a2 = z_pol_ring(neg_pol_a2_list)
            pol_g2 = pol_a2*self.h #Es calcula la convolucio.
            pol_g2 = self.r_ring_mod_q(pol_g2).lift().change_ring(z_pol_ring) #Es prenen moduls p i (x^N-1)
                                                                              #i es centre liften els coeficients.
            g2_list = list(pol_g2)

            if len(g2_list) < self.N: #Si la llista dels seus coeficients es menor que N, s'afegeixen els 0 que calgui al
                                      #final.
                zeros = [0]*(self.N - len(g2_list))
                g2_list.extend(zeros)

            g2_k_coef_list = g2_list[:k] #Es prenen els primers k coeficients del polinomi g2.
            self.classify(pol_a2_list, g2_k_coef_list, bins) #Es classifica el polinomi g2.

        g_list = [0]*self.N #S'inicialitzen a 0 l'index i el polinomi g.
        index = 0

        trobat = False

        while not trobat: #S'itera fins a trobar la parella de claus privades f i g.

            if len(bins[index])>1: #Es comproven tots els calaixos que disposin de mes d'un polinomi.

                pairs = itertools.combinations(bins[index],2) #Es calculen totes les possibles combinacions dels polinomis
                                                              #d'un mateix calaix.

                i=0
                pairs_list =list(pairs)

                while not trobat and i<len(pairs_list): #S'itera fins a trobar la parella de claus valides o esgotar totes
                                                        #les combinacions d'un calaix.
                    f_list = pairs_list[i][0] + pairs_list[i][1] #Es genera el candidat a clau privada f.
                    pol_f = z_pol_ring(f_list)
                    pol_g = pol_f*self.h #Es calcula la convolucio.
                    pol_g = center_lift(self.r_ring_mod_q(pol_g)) #Es prenen moduls p i (x^N-1)
                                                                  #i es centre liften els coeficients.
                    g_list = list(pol_g)

                    if len(g_list) < self.N: #Si la llista dels seus coeficients es menor que N, s'afegeixen els 0 que calgui al
                                             #final.
                        zeros = [0]*(self.N - len(g_list))
                        g_list.extend(zeros)

                    trobat = is_ternary(g_list, self.N, self.d) #Es comprova si el polinomi g compleix les condicions de
                                                                #la clau privada.
                    i+=1

            index+=1

        self.f = self.r_ring(f_list) #S'emmagatemen les claus privades f i g i els inversos modul p i q de f.
        self.g = self.r_ring(g_list)

        self.f_inv_mod_q = inv_pol_mod_power_2(self.r_ring_mod_q(f_list), self.r_ring_mod_2(f_list), self.N, self.q)

        if not self.p.is_power_of(2):
            self.f_inv_mod_p = inv_pol_mod_prime(self.r_ring_mod_p(f_list))

        else:
            self.f_inv_mod_p = inv_pol_mod_power_2(self.r_ring_mod_p(f_list), self.r_ring_mod_2(f_list), self.N, self.p)



    def classify(self, pol_a1, g1_k_coef_list, bins): #Es classifica el polinomi entrat als calaixos segons els bits
                                                      #dels seus primers k coeficients.

        first_bit_list = [] #Llista que conte els primers bits dels k primers coeficients.
        first_modified_bit_list = [] #Llista que conte els primers bits dels k primers coeficients modificats.
        number_1 = 0
        number_2 = 0

        num_of_bits = ceil(log(self.q,2)) #Es calcula el nombre de bits que tindra el valor resultant.

        first_bit_list, first_modified_bit_list = self.convert_to_bits(g1_k_coef_list, num_of_bits) #Es calculen els
                                                                                                    #primers bits de cada
                                                                                                    #coeficient.

        for bit_1 in first_bit_list: #S'obte el nombre corresponent a concatenar els bits.
            number_1 = (number_1 << 1) | bit_1

        for bit_2 in first_modified_bit_list: #S'obte el nombre corresponent a concatenar els bits modificats.
            number_2 = (number_2 << 1) | bit_2

        bins[number_1].append(pol_a1) #S'afegeix el polinomi als calaixos.
        bins[number_2].append(pol_a1)



    def convert_to_bits(self, g1_k_coef_list, num_of_bits): #Es calculen els primers bits de cada coeficient.

        bits = [] #Llista que conte l'equivalent en binari d'un coeficient donat.
        modified_bits = [] #Llista que conte l'equivalent en binari d'un coeficient modificat donat.
        first_bit_list = [] #Llista que conte tots el primer bit de cada coeficient.
        first_modified_bit_list = [] #Llista que conte tots el primer bit modificat de cada coeficient.

        for coef in g1_k_coef_list: #S'itera per tots els coeficients.
            coef_int=Integer(coef)
            bits = [int(x) for x in '{:08b}'.format(coef_int)] #Es tradueix a binari el coeficient en questio.
            bits = bits[(8-num_of_bits):] #S'eliminen els 0 de l'inici.
            first_bit_list.append(bits[0]) #S'afegeix a la llista el bit mes representatiu.

            modified_coef = coef_int+1 #Es modifica el coeficient.
            modified_bits = [int(x) for x in '{:08b}'.format(modified_coef)]#Es tradueix a binari el coeficient en questio.
            modified_bits = bits[(8-num_of_bits):] #S'eliminen els 0 de l'inici.
            first_modified_bit_list.append(modified_bits[0]) #S'afegeix a la llista el bit mes representatiu.

        return first_bit_list, first_modified_bit_list
