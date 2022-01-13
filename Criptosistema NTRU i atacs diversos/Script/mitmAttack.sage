load("auxMethods.sage")

import logging
import random
import numpy as np
import itertools
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
    bins_a1= None
    bins_a2 = None
    trobat = False
    g_list = None
    f_list = None
    k = None



    def __init__(self, *args): #Inicialitzacio dels parametres i dels anells que es requereixen per fer funcionar l'atac al KRP basat en reticles sense llegir un fitxer.

        if len(args)==1:
            f = open(args[0], "r") #Es llegeix el fitxer que conte els parametres publics.
            params = f.read()
            f.close()

            param_list = [] #Llista que conte els parametres publics en format string.
            param_list  = params.split(";")
            param_list[-1] = param_list[-1].split(",") #Es separen els coeficients de la clau publica h llegida del fitxer.
            h_int_list = [int(c) for c in param_list[-1]] #Es passen de coeficients de tipus string a tipus int.

            z_pol_ring.<x> = PolynomialRing(ZZ)
            self.N = Integer(param_list[0]) #S'estableixen els nous valors dels parametres publics.
            self.p = Integer(param_list[1])
            self.q = Integer(param_list[2])
            self.d = Integer(param_list[3])
            self.h = z_pol_ring(h_int_list)

        else:
            try:
                assert (args[0] - 1)/2 >= args[3] #Comprovem la condició 1 dels parametres.

            except AssertionError:
                print('Cal que (N - 1)/2 >= d')
                sys.exit(1)

            try:
                assert args[2] > (6*args[3] + 1)*args[1] #Comprovem la condició 2 dels parametres.

            except AssertionError:
                print('Cal que q > (6*d + 1)*p')
                sys.exit(1)

            try:
                assert args[2].is_power_of(2) #Comprovem que q sigui potència de 2.

            except AssertionError:
                print('q ha de ser potència de 2.')
                sys.exit(1)

            self.N = args[0]
            self.p = args[1]
            self.q = args[2]
            self.d = args[3]

        aux_ring.<x> = PolynomialRing(ZZ)
        self.r_ring = aux_ring.quotient(x^self.N - 1)

        aux_ring.<x> = PolynomialRing(ZZ.quotient(self.p))
        self.r_ring_mod_p = aux_ring.quotient(x^self.N - 1)

        aux_ring.<x> = PolynomialRing(ZZ.quotient(self.q))
        self.r_ring_mod_q = aux_ring.quotient(x^self.N - 1)

        aux_ring.<x> = PolynomialRing(ZZ.quotient(2))
        self.r_ring_mod_2 = aux_ring.quotient(x^self.N - 1)



    def decryption(self, e): #Es desencripta el missatge entrat per parametre.

        z_pol_ring.<x> = PolynomialRing(ZZ)
        a = z_pol_ring(list(self.f))*e #Es calcula la primera de les dues convolucions.
        a = center_lift(self.r_ring_mod_q(a)) #Es prenen moduls q i (x^N-1) i es centre liften els coeficients.
        b = z_pol_ring(list(a))*self.f_inv_mod_p #Es calcula la segona de les dues convolucions.
        b = center_lift(self.r_ring_mod_p(b)) #Es prenen moduls p i (x^N-1) i es centre liften els coeficients.

        return b



    def decryption_no_file(self, e): #Desencripta missatges sense lectura de fitxers.

        b = self.decryption(e)
        b_list = list(b)

        if (len(b_list)) < self.N: #Si el tamany de la llista es menor que N, cal afegir els 0 restants al final.
            zeros = [0]*(self.N - len(b_list))
            b_list.extend(zeros)

        return b_list



    def decrypt_message(self, list_of_encrypted): #Es desencripten el conjunt de submissatges encriptats.

        decrypted_list = []

        for i in range(len(list(list_of_encrypted))): #Per cadascun dels submissatges es crida a la funcio que els desencripta.
            decrypted_list.append(list(self.decryption(list_of_encrypted[i])))

        return decrypted_list #Es retorna una llista amb tots els submissatges desencriptats.



    def decrypt_file(self, e_filename, d_filename): #Donats dos fitxers, llegeix els coeficients dels polinomis corresponents al text encriptat i escriu el missatge desencriptat al segon fitxer.

        decrypted_lines_list = [] #Llista on els elements son les linies de text desencriptades.
        b_list = [] #Llista que conte subllistes corresponents als bits dels submissatges d'una linia de text desencriptada.
        b_def_list = [] #Llista que conte els bits dels submissatges unificats d'una linia de text desencriptada.
        pol_lines_list = [] #Llista que conte subllistes referents als coeficients dels missatges encriptats per cada linia de text.
        pol_lines_list = read_encrypted_file(e_filename, self.N) #Es llegeix el text encriptat del fitxer corresponent.

        for line in pol_lines_list: #Es desencripten les linies del fitxer una a una.
            b_list = self.decrypt_message(line) #Es desencripta una linia.
            add_zeros(b_list, self.N) #S'afegeixen els 0 que calgui al final de cada polinomi.
            b_def_list = [bit for polynomial in b_list for bit in polynomial] #S'unifiquen els coeficients dels polinomis a una sola linia.
            dec_text = "".join(map(chr, np.packbits(np.array(list(b_def_list))))) #Es transformen els bits al text en string corresponent.
            decrypted_lines_list.append(dec_text+"\n") #S'afegeix un salt de linia al final de cada linia.

        f = open(d_filename, "w") #S'escriu el text desencriptat al fitxer corresponent.
        f.writelines(decrypted_lines_list)
        f.close()



    def attack_file(self, e_filename, a_filename): #S'ataca el primer fitxer entrat per parametre i es guarda el missatge al segon fitxer entrat per parametre.

        try:
            self.collision(self.h)
            self.obtain_priv_keys() #S'obtenen les claus privades mitjançant l'atac.
            self.decrypt_file(e_filename, a_filename) #Es desencripta el fitxer.

        except (ZeroDivisionError, TypeError):
            print("Impossible obtenir les claus privades.")

        except MemoryError:
            print("Espai de memòria insuficient: impossible obtenir les claus privades.")



    def collision(self, h): #S'obtenen les claus privades seguint l'atac per trobada a mig cami.

        self.h = h

        try:
            self.g_list =[0]*self.N #Llista que conte els coeficients de la clau privada g.
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
            a1_iter = sympy.utilities.iterables.multiset_permutations(a1_list) #S'obtenen totes les possibles permutacions de la primera part de la clau privada f.
            a2_iter = sympy.utilities.iterables.multiset_permutations(a2_list) #S'obtenen totes les possibles permutacions de la segona part de la clau privada f.
            self.k = ceil(log((factorial(self.N)/(factorial(self.d/2+1)*factorial(self.d/2)*factorial(self.N-self.d-1))))) #Es calcula el nombre de calaixos que s'usaran per classificar.
            self.bins_a1 = [[] for Null in range(pow(2,self.k))] #Llista de tots els calaixos.
            self.bins_a2 = [[] for Null in range(pow(2,self.k))] #Llista de tots els calaixos.

            for a1_p in a1_iter: #Es classifiquen totes les possibles primeres parts del polinomi f
                addzeros_1 = [0]*(self.N-a1_N)
                pol_a1_list = addzeros_1 + list(a1_p)
                pol_a1 = z_pol_ring(pol_a1_list)
                pol_g1 = pol_a1*self.h  #Es calcula la convolucio.
                pol_g1 = self.r_ring_mod_q(pol_g1).lift().change_ring(z_pol_ring) #Es prenen moduls p i (x^N-1) i es centre liften els coeficients.
                g1_list = list(pol_g1)

                if len(g1_list) < self.N: #Si la llista dels seus coeficients es menor que N, s'afegeixen els 0 que calgui al final.
                    zeros = [0]*(self.N - len(g1_list))
                    g1_list.extend(zeros)

                g1_k_coef_list = g1_list[:self.k] #Es prenen els primers k coeficients del polinomi g1.
                self.classify(pol_a1_list, g1_k_coef_list, True) #Es classifica el polinomi g1.

            for a2_p in itertools.takewhile(self.checking, a2_iter): #Es classifiquen totes les possibles segones parts del polinomi f
                pass

        except MemoryError:
            raise



    def checking(self,a2_p):

        try:
            a2_N = floor(self.N/2)
            index = 0
            self.trobat = False
            z_pol_ring.<x> = PolynomialRing(ZZ)
            addzeros_2 = [0]*(self.N-a2_N)
            pol_a2_list = list(a2_p) + addzeros_2
            neg_pol_a2_list = [ x for x in pol_a2_list] #Es converteix el polinomi en negatiu.
            pol_a2 = z_pol_ring(neg_pol_a2_list)
            pol_g2 = -pol_a2*self.h #Es calcula la convolucio.
            pol_g2 = self.r_ring_mod_q(pol_g2).lift().change_ring(z_pol_ring) #Es prenen moduls p i (x^N-1) i es centre liften els coeficients.
            g2_list = list(pol_g2)

            if len(g2_list) < self.N: #Si la llista dels seus coeficients es menor que N, s'afegeixen els 0 que calgui al final.
                zeros = [0]*(self.N - len(g2_list))
                g2_list.extend(zeros)

            g2_k_coef_list = g2_list[:self.k] #Es prenen els primers k coeficients del polinomi g2.
            bins_list = []
            bins_list = self.classify(pol_a2_list, g2_k_coef_list, False) #Es classifica el polinomi g2.

            while not self.trobat and index <len(bins_list): #S'itera fins a trobar la parella de claus privades f i g.

                if len(self.bins_a1[bins_list[index]])>0: #Es comproven tots els calaixos que disposin de mes d'un polinomi.
                    i2=0

                    while not self.trobat and i2 < len(self.bins_a1[bins_list[index]]):
                        pol_f_1 = z_pol_ring(self.bins_a1[bins_list[index]][i2])
                        pol_f_2 = z_pol_ring(pol_a2_list)
                        pol_f = pol_f_1+pol_f_2
                        self.f_list = list(pol_f)
                        pol_g = pol_f*self.h #Es calcula la convolucio.
                        pol_g = center_lift(self.r_ring_mod_q(pol_g)) #Es prenen moduls p i (x^N-1) i es centre liften els coeficients.
                        self.g_list = list(pol_g)

                        if len(self.g_list) < self.N: #Si la llista dels seus coeficients es menor que N, s'afegeixen els 0 que calgui al final.
                            zeros = [0]*(self.N - len(self.g_list))
                            self.g_list.extend(zeros)

                        self.trobat = is_ternary(self.g_list, self.N, self.d) #Es comprova si el polinomi g compleix les condicions de la clau privada.
                        i2+=1

                index+=1

            return not self.trobat

        except MemoryError:
            raise



    def obtain_priv_keys(self):

        z_pol_ring.<x> = PolynomialRing(ZZ)

        if self.trobat:
            try:
                self.f = self.r_ring(self.f_list) #S'emmagatemen les claus privades f i g i els inversos modul p i q de f.
                self.g = self.r_ring(self.g_list)
                self.f_inv_mod_q = inv_pol_mod_power_2(self.r_ring_mod_q(self.f_list), self.r_ring_mod_2(self.f_list), self.N, self.q)

                if not self.p.is_power_of(2):
                    self.f_inv_mod_p = inv_pol_mod_prime(self.r_ring_mod_p(self.f_list))

                else:
                    self.f_inv_mod_p = inv_pol_mod_power_2(self.r_ring_mod_p(self.f_list), self.r_ring_mod_2(self.f_list), self.N, self.p)

            except ZeroDivisionError:
                raise

            except MemoryError:
                raise



    def classify(self, pol_a1, g1_k_coef_list, b): #Es classifica el polinomi entrat als calaixos segons els bits dels seus primers k coeficients.

        try:
            first_bit_list = [] #Llista que conte els primers bits dels k primers coeficients.
            first_modified_bit_list = [] #Llista que conte els primers bits dels k primers coeficients modificats.
            number_1 = 0
            number_2 = 0
            num_of_bits = ceil(log(self.q,2)) #Es calcula el nombre de bits que tindra el valor resultant.
            first_bit_list, first_modified_bit_list, first_modified_bit_list_2 = self.convert_to_bits(g1_k_coef_list, num_of_bits) #Es calculen els primers bits de cada coeficient.
            list_of_first_bit_list = []
            list_of_first_bit_list.append(first_bit_list)

            for i in range(0,len(first_bit_list)):

                if first_bit_list[i] != first_modified_bit_list[i]:
                    l = len(list_of_first_bit_list)

                    for e in range(0,l):
                        new_list = list_of_first_bit_list[e].copy()
                        new_list[i] = first_modified_bit_list[i]
                        list_of_first_bit_list.append(new_list)


                if first_bit_list[i] != first_modified_bit_list_2[i]:
                    l = len(list_of_first_bit_list)

                    for e in range(0,l):
                        new_list = list_of_first_bit_list[e].copy()
                        new_list[i] = first_modified_bit_list_2[i]
                        list_of_first_bit_list.append(new_list)

            counter = 0
            bins_list = []

            for list_of_bits in list_of_first_bit_list:
                counter+=1

                for bit_1 in list_of_bits: #S'obte el nombre corresponent a concatenar els bits.
                    number_1 = (number_1 << 1) | bit_1

                if b:
                    self.bins_a1[number_1].append(pol_a1) #S'afegeix el polinomi als calaixos.

                else:
                    bins_list.append(number_1)
                    self.bins_a2[number_1].append(pol_a1) #S'afegeix el polinomi als calaixos.

                number_1 = 0

            if not b:
                return bins_list

        except MemoryError:
            raise



    def convert_to_bits(self, g1_k_coef_list, num_of_bits): #Es calculen els primers bits de cada coeficient.

        bits = [] #Llista que conte l'equivalent en binari d'un coeficient donat.
        modified_bits = [] #Llista que conte l'equivalent en binari d'un coeficient modificat donat.
        modified_bits_2 = [] #Llista que conte l'equivalent en binari d'un coeficient modificat donat.
        first_bit_list = [] #Llista que conte tots el primer bit de cada coeficient.
        first_modified_bit_list = [] #Llista que conte tots el primer bit modificat de cada coeficient.
        first_modified_bit_list_2 = [] #Llista que conte tots el primer bit modificat de cada coeficient.

        for coef in g1_k_coef_list: #S'itera per tots els coeficients.
            coef_int=Integer(coef)
            bits = [int(x) for x in '{:020b}'.format(coef_int)] #Es tradueix a binari el coeficient en questio.
            bits = bits[(20-num_of_bits):] #S'eliminen els 0 de l'inici.
            first_bit_list.append(bits[0]) #S'afegeix a la llista el bit mes representatiu.
            modified_coef = (coef_int+1)%self.q #Es modifica el coeficient.
            modified_bits = [int(x) for x in '{:020b}'.format(modified_coef)]#Es tradueix a binari el coeficient en questio.
            modified_bits = modified_bits[(20-num_of_bits):] #S'eliminen els 0 de l'inici.
            first_modified_bit_list.append(modified_bits[0]) #S'afegeix a la llista el bit mes representatiu.
            modified_coef_2 = (coef_int-1)%self.q #Es modifica el coeficient.
            modified_bits_2 = [int(x) for x in '{:020b}'.format(modified_coef_2)]#Es tradueix a binari el coeficient en questio.
            modified_bits_2 = modified_bits_2[(20-num_of_bits):] #S'eliminen els 0 de l'inici.
            first_modified_bit_list_2.append(modified_bits_2[0]) #S'afegeix a la llista el bit mes representatiu.

        return first_bit_list, first_modified_bit_list, first_modified_bit_list_2
