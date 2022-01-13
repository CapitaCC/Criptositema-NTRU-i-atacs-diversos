load("auxMethods.sage")

import logging
import random
import numpy as np
import itertools



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



    def attack_file(self, e_filename, d_filename): #S'ataca el primer fitxer entrat per parametre i es guarda el missatge al segon fitxer entrat per parametre

        decrypted_lines_list = [] #Llista on els elements son les linies de text desencriptades.
        a_list = [] #Llista que conte subllistes corresponents als bits dels submissatges d'una linia de text desencriptada.
        a_def_list = [] #Llista que conte els bits dels submissatges unificats d'una linia de text desencriptada.
        pol_lines_list = [] #llista que conte subllistes referents als coeficients dels missatges encriptats per cada linia de text.
        gc = gcd(self.q,self.p)
        exit = False

        if gc == 1:
            print("{} i {} coprimers, atac impossible de realitzar.".format(self.p, self.q))

        else:
            pol_lines_list = read_encrypted_file(e_filename, self.N) #Es llegeix el text encriptat del fitxer corresponent.

            for line in pol_lines_list: #Es desencripten les linies del fitxer una a una.
                a_list = self.attack_message(line, gc) #Es desencripta una linia.
                add_zeros(a_list, self.N) #S'afegeixen els 0 que calgui al final de cada polinomi.
                a_def_list = [bit for polynomial in a_list for bit in polynomial] #S'unifiquen els coeficients dels polinomis a una sola linia.
                dec_text = "".join(map(chr, np.packbits(np.array(list(a_def_list))))) #Es transformen els bits al text en string corresponent.
                decrypted_lines_list.append(dec_text+"\n") #S'afegeix un salt de linia al final de cada linia.

            exit = True

        f = open(d_filename, "w") #S'escriu el text desencriptat al fitxer corresponent.
        f.writelines(decrypted_lines_list)
        f.close()

        return exit


    def attack_message(self, list_of_encrypted, gc): #S'ataquen el conjunt de submissatges encriptats.

        attacked_list = []

        for i in range(len(list(list_of_encrypted))): #Per cadascun dels submissatges es crida a la funcio que els ataca.
            attacked_list.append(list(self.attack(list_of_encrypted[i], gc)))

        return attacked_list #Es retorna una llista amb tots els submissatges atacats.



    def attack(self, e, gc): #S'ataca el missatge entrat per parametre.

        z_pol_ring.<x> = PolynomialRing(ZZ)
        aux_ring.<x> = PolynomialRing(ZZ.quotient(Integer(gc)))
        r_ring_mod_c = aux_ring.quotient(x^self.N - 1)
        b = center_lift(r_ring_mod_c(list(e)))

        return b



    def attack_no_file(self, e): #S'ataca el missatge entrat per parametre sense llegir de fitxers.

        gc = gcd(self.q, self.p)

        if gc == 1:
            print("{} i {} coprimers, atac impossible de realitzar.".format(self.p, self.q))
            b = [0]

        else:
            z_pol_ring.<x> = PolynomialRing(ZZ)
            aux_ring.<x> = PolynomialRing(ZZ.quotient(Integer(gc)))
            r_ring_mod_c = aux_ring.quotient(x^self.N - 1)
            b = list(center_lift(r_ring_mod_c(list(e))))

            if (len(b)) < self.N: #Si el tamany de la llista es menor que N, cal afegir els 0 restants al final.
                zeros = [0]*(self.N - len(b))
                b.extend(zeros)

        return b