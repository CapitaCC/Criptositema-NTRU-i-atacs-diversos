load("auxMethods.sage")

import logging
import random
import numpy as np
import itertools
import time
import base64
import sympy



class Param_vulnerability:

    N = None
    p = None
    q = None
    d = None
    r = None
    r_ring = None
    r_ring_mod_p = None
    r_ring_mod_q = None
    r_ring_mod_2 = None




    def __init__(self, filename): #Inicialitzacio dels parametres i dels anells que es requereixen per fer funcionar
                                    #l'atac a la vulnerabilitat dels parametres.

        f = open(filename, "r") #Es llegeix el fitxer que conte els parametres publics.
        params = f.read()
        f.close()

        param_list = [] #Llista que conte els parametres publics en format string.

        param_list  = params.split("A")

        self.N = Integer(param_list[0]) #S'estableixen els nous valors dels parametres publics.
        self.p = Integer(param_list[1])
        self.q = Integer(param_list[2])
        self.d = Integer(param_list[3])
        aux_ring.<x> = PolynomialRing(ZZ)
        self.r_ring = aux_ring.quotient(x^self.N - 1)

        aux_ring.<x> = PolynomialRing(ZZ.quotient(self.p)) #PolynomialRing(QuotientRing(ZZ, ZZ.ideal(self.p)))
        self.r_ring_mod_p = aux_ring.quotient(x^self.N - 1) #QuotientRing(aux_ring, aux_ring.ideal(x^self.N - 1))

        aux_ring.<x> = PolynomialRing(ZZ.quotient(self.q)) #PolynomialRing(QuotientRing(ZZ, ZZ.ideal(self.q)))
        self.r_ring_mod_q = aux_ring.quotient(x^self.N - 1) #QuotientRing(aux_ring, aux_ring.ideal(x^self.N - 1))

        aux_ring.<x> = PolynomialRing(ZZ.quotient(2))

    def attack_file(self, e_filename, d_filename):

        decrypted_lines_list = [] #Llista on els elements son les linies de text desencriptades.
        a_list = [] #Llista que conte subllistes corresponents als bits dels submissatges d'una linia de text
                    #desencriptada.
        a_def_list = [] #Llista que conte els bits dels submissatges unificats d'una linia de text desencriptada.
        pol_lines_list = [] #llista que conte subllistes referents als coeficients dels missatges encriptats per cada
                            #linia de text.

        gc = gcd(self.q,self.p)

        if gc == 1:

            print("{} i {} coprimers, atac impossible de realitzar.".format(self.p, self.q))
            return False

        else:

            pol_lines_list = read_encrypted_file(e_filename, self.N) #Es llegeix el text encriptat del fitxer corresponent.

            for line in pol_lines_list: #Es desencripten les linies del fitxer una a una.
                a_list = self.attack_message(line, gc) #Es desencripta una linia.
                add_zeros(a_list, self.N) #S'afegeixen els 0 que calgui al final de cada polinomi.
                a_def_list = [bit for polynomial in a_list for bit in polynomial] #S'unifiquen els coeficients dels polinomis
                                                                                  #a una sola linia.
                dec_text = "".join(map(chr, np.packbits(np.array(list(a_def_list))))) #Es transformen els bits al text en
                                                                                      #string corresponent.
                decrypted_lines_list.append(dec_text+"\n") #S'afegeix un salt de linia al final de cada linia.

        f = open(d_filename, "w") #S'escriu el text desencriptat al fitxer corresponent.
        f.writelines(decrypted_lines_list)
        f.close()
        return True


    def attack_message(self, list_of_encrypted, gc): #S'ataquen el conjunt de submissatges encriptats.

        attacked_list = []

        for i in range(len(list(list_of_encrypted))): #Per cadascun dels submissatges es crida a la funcio que els
                                                      #ataca.

            attacked_list.append(list(self.attack(list_of_encrypted[i], gc)))

        return attacked_list #Es retorna una llista amb tots els submissatges atacats.


    def attack(self, e, gc): #S'ataca el missatge entrat per parametre.

        z_pol_ring.<x> = PolynomialRing(ZZ)
        aux_ring.<x> = PolynomialRing(ZZ.quotient(Integer(gc))) #PolynomialRing(QuotientRing(ZZ, ZZ.ideal(self.q)))
        r_ring_mod_c = aux_ring.quotient(x^self.N - 1) #QuotientRing(aux_ring, aux_ring.ideal(x^self.N - 1))


        b = center_lift(r_ring_mod_c(list(e)))

        return b