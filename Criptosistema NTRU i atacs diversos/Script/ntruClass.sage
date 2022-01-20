load("auxMethods.sage")

import logging
import random
import numpy as np
import itertools



class Ntru:

    N = None
    p = None
    q = None
    d = None
    f = None
    g = None
    h = None
    r = None
    f_inv_mod_p = None
    f_inv_mod_q = None
    r_ring = None
    r_ring_mod_p = None
    r_ring_mod_q = None
    r_ring_mod_2 = None



    def __init__(self, N, p, q, d, info): #Inicialització dels parametres i dels anells que es requereixen per fer funcionar el criptosistema NTRU.

        try:
            assert (N - 1)/2 >= d #Comprovem la condició 1 dels parametres.

        except AssertionError:
            print('Cal que (N - 1)/2 >= d')
            info = False

        try:
            assert q > (6*d + 1)*p #Comprovem la condició 2 dels parametres.

        except AssertionError:
            print('Cal que q > (6*d + 1)*p')
            info = False

        try:
            assert q.is_power_of(2) #Comprovem que q sigui potència de 2.

        except AssertionError:
            print('q ha de ser potència de 2.')
            info = False

        self.N = N
        self.p = p
        self.q = q
        self.d = d

        aux_ring.<x> = PolynomialRing(ZZ)
        self.r_ring = aux_ring.quotient(x^self.N - 1)

        aux_ring.<x> = PolynomialRing(ZZ.quotient(self.p))
        self.r_ring_mod_p = aux_ring.quotient(x^self.N - 1)

        aux_ring.<x> = PolynomialRing(ZZ.quotient(self.q))
        self.r_ring_mod_q = aux_ring.quotient(x^self.N - 1)

        aux_ring.<x> = PolynomialRing(ZZ.quotient(2))
        self.r_ring_mod_2 = aux_ring.quotient(x^self.N - 1)

        if info:
            print("Criptosistema NTRU instanciat amb els paràmetres N={}, p={}, q={}, d={}".format(N, p, q, d))



    def gen_priv_keys(self): #Es calculen les claus privades f i g, aixi com tambe els polinomis f_inv_mod_p i f_inv_mod_q.

        trobat = False
        g_list = gen_ternary_pol(False, self.N, self.d) #S'obtenen els coeficients del polinomi g.
        self.g = self.r_ring(g_list)

        while not trobat: #Mentre no s'hagi trobat un f que disposi d'inversos modul p i q, se'n busquen de nous.
            f_list = gen_ternary_pol(True, self.N, self.d) #S'obtenen els coeficients del polinomi g.
            self.f = self.r_ring(f_list)

            if not self.p.is_power_of(2): #Cas en que p no es potencia de 2.
                aux_f_p = self.r_ring_mod_p(f_list)

                if aux_f_p.is_unit(): #En cas que f tingui invers modul p, es calcula.
                    self.f_inv_mod_p = inv_pol_mod_prime(aux_f_p)
                    aux_f_q = self.r_ring_mod_2(f_list)

                    if aux_f_q.is_unit(): #En cas que f tingui invers modul q, s'aplica l'algorisme iteratiu de Newton per obtenir l'invers modul q de f.
                        self.f_inv_mod_q = inv_pol_mod_power_2(self.r_ring_mod_q(f_list), aux_f_q, self.N, self.q)
                        trobat = True #S'indica que s'ha trobat una clau privada f valida i els seus inversos moduls p i q.

            else: #Cas en que p es potencia de 2.
                aux_f_p = self.r_ring_mod_2(f_list)

                if aux_f_p.is_unit(): #En cas que f tingui invers modul p, es calcula.
                    self.f_inv_mod_p = inv_pol_mod_power_2(self.r_ring_mod_p(f_list), aux_f_p, self.N, self.p)
                    aux_f_q = self.r_ring_mod_2(f_list)

                    if aux_f_q.is_unit(): #En cas que f tingui invers modul q, s'aplica l'algorisme iteratiu de Newton per obtenir l'invers modul q de f.
                        self.f_inv_mod_q = inv_pol_mod_power_2(self.r_ring_mod_q(f_list), aux_f_q, self.N, self.q)
                        trobat = True #S'indica que s'ha trobat una clau priva f valida i els seus inversos moduls p i q.



    def gen_pub_key(self): #Es calcula la clau publica h.

        z_pol_ring.<x> = PolynomialRing(ZZ)
        h_z = z_pol_ring(list(self.f_inv_mod_q))*z_pol_ring(list(self.g)) #El producte entre f i g es fa a l'anell dels polinomis amb coeficients enters.
        self.h = center_lift(self.r_ring_mod_q(h_z)) #Es prenen moduls q i (x^N-1) i es centre liften els coeficients.



    def gen_keys(self, pub_filename, priv_filename): #Es calculen tant les claus privades com la publica.

        self.gen_priv_keys()
        self.gen_pub_key()
        store_param_and_pub_key(pub_filename, self.N, self.p, self.q, self.d, self.h)
        store_param_and_priv_keys(priv_filename, self.N, self.p, self.q, self.d, self.f, self.g)



    def set_param_and_pub_key(self, filename): #Estableix els parametres i la clau publica segons els valors del fitxer indicat.

        f = open(filename, "r") #Es llegeix el fitxer que conte els parametres i la clau publica.
        params = f.read()
        f.close()

        param_list = [] #Llista que conte els parametres publics en format string.
        h_int_list = [] #Llista que conte els coeficients de la clau publica h en format int.
        param_list  = params.split(";") #Es separen els parametres publics llegits del fitxer de text.
        param_list[-1] = param_list[-1].split(",") #Es separen els coeficients de la clau publica h llegida del fitxer.
        h_int_list = [int(c) for c in param_list[-1]] #Es passen de coeficients de tipus string a tipus int.

        z_pol_ring.<x> = PolynomialRing(ZZ)
        self.N = Integer(param_list[0]) #S'estableixen els nous valors dels parametres publics.
        self.p = Integer(param_list[1])
        self.q = Integer(param_list[2])
        self.d = Integer(param_list[3])

        aux_ring.<x> = PolynomialRing(ZZ)
        self.r_ring = aux_ring.quotient(x^self.N - 1)

        aux_ring.<x> = PolynomialRing(ZZ.quotient(self.p))
        self.r_ring_mod_p = aux_ring.quotient(x^self.N - 1)

        aux_ring.<x> = PolynomialRing(ZZ.quotient(self.q))
        self.r_ring_mod_q = aux_ring.quotient(x^self.N - 1)

        aux_ring.<x> = PolynomialRing(ZZ.quotient(2))
        self.r_ring_mod_2 = aux_ring.quotient(x^self.N - 1)
        self.h = z_pol_ring(h_int_list)



    def set_param_and_priv_keys(self, filename): #Estableix els parametres i les claus privades segons els valors del fitxer indicat.

        f = open(filename, "r") #Es llegeix el fitxer que conte els parametres i la clau publica.
        params = f.read()
        f.close()

        param_list = [] #Llista que conte els parametres publics en format string.
        f_int_list = [] #Llista que conte els coeficients de la clau privada f en format int.
        g_int_list = [] #Llista que conte els coeficients de la clau privada g en format int.
        param_list  = params.split(";") #Es separen els parametres publics llegits del fitxer de text.
        param_list[-2] = param_list[-2].split(",") #Es separen els coeficients de la clau privada f llegida del fitxer.
        param_list[-1] = param_list[-1].split(",") #Es separen els coeficients de la clau privada g llegida del fitxer.
        f_int_list = [int(c) for c in param_list[-2]] #Es passen de coeficients de tipus string a tipus int.
        g_int_list = [int(c) for c in param_list[-1]] #Es passen de coeficients de tipus string a tipus int.

        z_pol_ring.<x> = PolynomialRing(ZZ)
        self.N = Integer(param_list[0]) #S'estableixen els nous valors dels parametres publics.
        self.p = Integer(param_list[1])
        self.q = Integer(param_list[2])
        self.d = Integer(param_list[3])
        aux_ring.<x> = PolynomialRing(ZZ)
        self.r_ring = aux_ring.quotient(x^self.N - 1)

        aux_ring.<x> = PolynomialRing(ZZ.quotient(self.p))
        self.r_ring_mod_p = aux_ring.quotient(x^self.N - 1)

        aux_ring.<x> = PolynomialRing(ZZ.quotient(self.q))
        self.r_ring_mod_q = aux_ring.quotient(x^self.N - 1)

        aux_ring.<x> = PolynomialRing(ZZ.quotient(2))
        self.r_ring_mod_2 = aux_ring.quotient(x^self.N - 1)

        self.f = self.r_ring(f_int_list)
        self.g = self.r_ring(g_int_list)
        self.f_inv_mod_q = inv_pol_mod_power_2(self.r_ring_mod_q(f_int_list), self.r_ring_mod_2(f_int_list), self.N, self.q)

        if not self.p.is_power_of(2):
            self.f_inv_mod_p = inv_pol_mod_prime(self.r_ring_mod_p(f_int_list))

        else:
            self.f_inv_mod_p = inv_pol_mod_power_2(self.r_ring_mod_p(f_int_list), self.r_ring_mod_2(f_int_list), self.N, self.p)



    def encryption(self, m): #S'encripta el missatge entrat per parametre.

        try:
            assert len(m) <= self.N #Es comprova que el missatge tingui el tamany permes.

        except AssertionError:
            print('Missatge massa llarg')

        z_pol_ring.<x> = PolynomialRing(ZZ)
        m_pol = z_pol_ring(m)
        self.r = z_pol_ring(gen_ternary_pol(False, self.N, self.d)) #Es calcula el polinomi aleatori r per aquest missatge en concret.
        e = self.p*self.h*self.r + m_pol
        e = center_lift(self.r_ring_mod_q(e))

        return e



    def decryption(self, e): #Es desencripta el missatge entrat per parametre.

        z_pol_ring.<x> = PolynomialRing(ZZ)
        a = z_pol_ring(list(self.f))*e #Es calcula la primera de les dues convolucions.
        a = center_lift(self.r_ring_mod_q(a)) #Es prenen moduls q i (x^N-1) i es centre liften els coeficients.
        b = z_pol_ring(list(a))*self.f_inv_mod_p #Es calcula la segona de les dues convolucions.
        b = center_lift(self.r_ring_mod_p(b)) #Es prenen moduls p i (x^N-1) i es centre liften els coeficients.

        return b



    def encrypt_message(self, list_of_messages): #S'encripten el conjunt de submissatges.

        encrypted_list = []

        for i in range(len(list_of_messages)): #Per cadascun dels submissatges es crida a la funcio que els encripta.
            encrypted_list.append(self.encryption(list(list_of_messages[i])))

        return encrypted_list #Es retorna una llista amb tots els submissatges encriptats.



    def decrypt_message(self, list_of_encrypted): #Es desencripten el conjunt de submissatges encriptats.

        decrypted_list = []

        for i in range(len(list(list_of_encrypted))): #Per cadascun dels submissatges es crida a la funcio que els desencripta.
            decrypted_list.append(list(self.decryption(list_of_encrypted[i])))

        return decrypted_list #Es retorna una llista amb tots els submissatges desencriptats.



    def encrypt_file(self, m_filename, e_filename): #Donats dos fitxers, llegeix el text del primer i n'escriu els coeficients dels polinomis corresponents a la seva encriptacio al segon fitxer.

        write_lines_list = [] #Llista on cada element es correspon amb el text encriptat d'una linia.
        coef_list = [] #Llista on els elements son els coeficients (enters) dels polinomis corresponents a l'encriptacio d'una lina de text.
        string_list = [] #Llista on els elements son els coeficients (strings) dels polinomis corresponents a l'encriptacio d'una lina de text.
        lines_list = [] #Llista on els elements son cadascuna de les linies de text del fitxer a encriptar.
        m_list = [] #llista que conte els submissatges per cada linia de text a encriptar

        f = open(m_filename, "rb") #Es llegeix el fitxer de text a encriptar.
        lines_list = f.read().splitlines()
        f.close()

        for line in lines_list: #S'encripten les linies del fitxer una a una.
            m_bits = np.unpackbits(np.frombuffer(line, dtype=np.uint8)) #Es passa el text al seu equivalent en bits.
            m_bits = np.trim_zeros(m_bits,'b') #S'eliminen els 0 sobrants del final.
            split_message(list(m_bits), m_list, self.N) #Es divideix cada linia en submissatges per poder ser encriptats.
            e_list = self.encrypt_message(m_list) #S'encripten els submissatges.
            m_list = [] #Es torna a buidar la llista de submissatges.
            coef_list = pol_to_coef(e_list, self.N) #Es passa dels polinomis a la llista dels seus coeficients.
            string_list = [str(c) for c in coef_list] #Es passa dels coeficients enters a els coeficients en string.
            c_string = ",".join(string_list) #S'uneixen els coeficients en un sol string separat per comes.
            c_string +=";" #Per indicar el fi de la linia s'afegeix un caracter ;.
            write_lines_list.append(c_string) #S'afegeix el string resultant a la llista de linies encriptades.

        f = open(e_filename, "w") #S'escriuen les linies encriptades al fitxer que ha de contenir el missatge encriptat.
        f.writelines(write_lines_list)
        f.close()



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



    def decryption_no_file(self, e): #Desencripta missatges sense lectura de fitxers.

        b = self.decryption(e)
        b_list = list(b)

        if (len(b_list)) < self.N: #Si el tamany de la llista es menor que N, cal afegir els 0 restants al final.
            zeros = [0]*(self.N - len(b_list))
            b_list.extend(zeros)

        return b_list



    def get_pub_key(self): #Retorna la clau publica.

        return self.h



    def set_pub_key(self, new_h): #Assigna la clau publica el polinomi per paràmetre.

        self.h = new_h
