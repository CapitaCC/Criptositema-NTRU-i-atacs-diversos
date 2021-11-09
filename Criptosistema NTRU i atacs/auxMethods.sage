import logging
import random
import numpy as np
import itertools
def center_lift(pol): #Centra els coeficients del polinomi entrat per parametre a [-q/2,q/2].

        return pol.lift().map_coefficients(lambda c: c.lift_centered(), ZZ)


def gen_ternary_pol(isf, N, d): #Es genera un polinomi ternari amb d (o d+1) coeficients iguals a 1 i d
                                          #coeficients iguals a -1.

    key = [0]*N #Inicialitzem tots els N coeficients del polinomi a 0.
    index_list = [i for i in range(0, N)] #Prenem un llistat de tots els coeficients que valen 0.

    for j in range(d): #Per cadascuna de les d iteracions escollim dos coeficients buits i els omplim amb 1 i -1.
        index_1 = sage.misc.prandom.choice(index_list)
        index_list.remove(index_1)

        index_2 = sage.misc.prandom.choice(index_list)
        index_list.remove(index_2)

        key[index_1] = 1
        key[index_2] = -1

    if isf: #En cas que es tracti del polinomi f, s'afegeix un altre coeficient igual a 1.
        index_1 = sage.misc.prandom.choice(index_list)
        index_list.remove(index_1)

        key[index_1] = 1

    return key



def inv_pol_mod_prime(pol): #Donat un polinomi pertanyent a un anell quocient amb coeficients a un cos,
                                  #retorna l'invers del polinomi a aquest anell prenent els coeficients a (-p/2,p/2).
                                  #El polinomi retornat viu a l'anell dels polinomis am coeficients als enters.

    return center_lift(pol^(-1))



#Metode iteratiu de Newton
def inv_pol_mod_power_2(pol, aux, N, q): #Donat un polinomi pertanyent a un anell quocient amb coeficients als enters
                                         #modul q una potencia de 2, retorna l'invers del polinomi a aquest anell
                                         #prenent els coeficients a (-q/2,q/2).
                                         #El polinomi retornat viu a l'anell dels polinomis am coeficients als enters.

    index = 2
    inv_pol = inv_pol_mod_prime(aux) #Es calcula l'invers modul 2 del polinomi original.
    pol2 = center_lift(pol) #Es centre lifta el polinomi original per tenir els coeficients a (-q/2,q/2).

    while index < q: #S'itera per cada potència de 2 fins arribar a q aplicant el mètode iteratiu de Newton.
        index = index**2
        inv_pol = (2*inv_pol - pol2*(inv_pol)^(2))
        aux_ring.<x> = PolynomialRing(ZZ.quotient(index))
        aux_ring = aux_ring.quotient(x^N - 1)
        inv_pol = center_lift(aux_ring(inv_pol)) #Es centre lifta el polinomi resultant per tenir els coeficients
                                                 #a (-q/2,q/2).
    return inv_pol





def add_zeros(list_of_lists, N): #Afegeix 0 als darrers coeficients del polinomi que valen 0, ja que SageMath
                                    #no els comptabilitza, per tal que totes les llistes tinguin mida N.

    for l in list_of_lists: #Per cadascun dels llistats de coeficients afegim els 0 al final que calgui.

        if (len(l)) < N: #Si el tamany de la llista es menor que N, cal afegir els 0 restants al final.

            zeros = [0]*(N - len(l))
            l.extend(zeros)



def split_message(message, list_of_messages, N): #Es divideix el missatge en missatges per poder ser encriptats.

    if len(message) <= N: #En cas que la mida del missatge sigui menor que el maxim establert, es pot encriptar.

        list_of_messages.append(message)

    else: #En cas contrari, cal separar-lo en submissatges de mida N.

        list_of_messages.append(message[:N])
        split_message(message[N:], list_of_messages, N)


def store_param_and_pub_key(filename, N, p, q, d, h): #Emmagatzema els parametres i la clau publica al fitxer
                                                            #indicat.


    h_coef_list = [] #Llista que conte els coeficients de la clau publica h.
    param_list = []  #Llista que conte els parametres publics en format string.

    h_coef_list = list(h) #S'obte la llista de coeficients de la clau publica h.

    if len(h_coef_list) < N: #En cas que la llista de coeficients sigui menor que N, cal afegir-hi 0 al final.
        zeros = [0]*(N - len(h_coef_list))
        h_coef_list.extend(zeros)

    string_list = [str(c) for c in h_coef_list] #Es passa dels coeficients enters a els coeficients en string.
    h_c_string = ",".join(string_list) #S'uneixen els coeficients en un sol string separat per comes.

    param_list.append(str(N)) #S'afegeixen els parametres publics en format string.
    param_list.append(str(p))
    param_list.append(str(q))
    param_list.append(str(d))
    param_list.append(h_c_string)

    param_text = "A".join(param_list) #S'uneixen els parametres publics en un sol string separat per la lletra "A".

    f = open(filename, "w") #S'escriu el text en el fitxer corresponent.
    f.write(param_text)
    f.close()


def store_param_and_priv_keys(filename, N, p, q, d, f, g): #Emmagatzema els parametres i les claus privades al
                                                                 #fitxer indicat.


    f_coef_list = [] #Llista que conte els coeficients de la clau privada f.
    g_coef_list = [] #Llista que conte els coeficients de la clau privada g.
    param_list = []  #Llista que conte els parametres publics en format string.

    f_coef_list = list(f) #S'obte la llista de coeficients de la clau privada f.
    g_coef_list = list(g) #S'obte la llista de coeficients de la clau privada g.

    if len(f_coef_list) < N: #En cas que la llista de coeficients sigui menor que N, cal afegir-hi 0 al final.
        zeros = [0]*(N - len(f_coef_list))
        f_coef_list.extend(zeros)

    if len(g_coef_list) < N: #En cas que la llista de coeficients sigui menor que N, cal afegir-hi 0 al final.
        zeros = [0]*(N - len(g_coef_list))
        g_coef_list.extend(zeros)

    string_list = [str(c) for c in f_coef_list] #Es passa dels coeficients enters a els coeficients en string.
    f_c_string = ",".join(string_list) #S'uneixen els coeficients en un sol string separat per comes.

    string_list = [str(c) for c in g_coef_list] #Es passa dels coeficients enters a els coeficients en string.
    g_c_string = ",".join(string_list) #S'uneixen els coeficients en un sol string separat per comes.

    param_list.append(str(N)) #S'afegeixen els parametres publics en format string.
    param_list.append(str(p))
    param_list.append(str(q))
    param_list.append(str(d))
    param_list.append(f_c_string)
    param_list.append(g_c_string)

    param_text = "A".join(param_list) #S'uneixen els parametres publics en un sol string separat per la lletra "A".

    f = open(filename, "w") #S'escriu el text en el fitxer corresponent.
    f.write(param_text)
    f.close()



def read_encrypted_file(filename, N): #Llegeix el fixer encriptat i en converteix el text a una llista que conte
                                         #subllistes referents als polinomis corresponents als missatges encriptats
                                         #per cada linia de text.

    f = open(filename, "r") #Es llegeix el fitxer que conte el missatge encriptat.
    coef = f.read()
    f.close()

    lines_list  = [] #Llista que conte el text corresponent a cadascuna de les linies encriptades.
    pol_lines_list  = [] #llista que conte subllistes referents als coeficients dels missatges encriptats per cada
                         #linia de text.
    int_coef_line_list = [] #Llista que conte els coeficients (enters) dels polinomis corresponents a una linia del
                            #text encriptat.
    pol_line_list  = [] #Llista que conte els polinomis corresponents al text encriptat d'una linia de text.

    lines_list  = coef.split("A") #Es separen els coeficients llegits del fitxer per línies de text.
    lines_list.pop() #S'elimina l'ultim element ja que es correspon amb una linia buida.

    z_pol_ring.<x> = PolynomialRing(ZZ) #Es crea l'anell de polinomis amb coeficients enters.

    for line in lines_list: #Per cada linia de text es transformen els coeficients en els polinomis corresponents.

        int_coef_line_list = [int(c) for c in line.split(',')] #Es passen de coeficients de tipus string a tipus int.

        i=0 #S'inicia el comptador a 0.

        for coef in int_coef_line_list: #S'itera pels coeficients de cada linia de text.

            if i%(N) == 0: #Cada N coeficiens cal crear un nou polinomi i afegir-lo a la llista.

                pol_line_list.append(z_pol_ring(int_coef_line_list[i:i+N]))
            i+=1

        pol_lines_list.append(pol_line_list) #S'afegeixen els polinomis corresponents a una linia de text.
        pol_line_list = [] #Es buida la llista corresponent als polinomis d'una linia de text.

    return pol_lines_list


def pol_to_coef(pol_list, N): #Es converteixen els polinomis en llistes que contenen els seus coeficients.

    coef_list = [] #Llista de subllistes que conte els coeficients ordenats dels polinomis.

    for pol in pol_list: #Es trasforma cada polinomi a una llista dels seus coeficients.

        c_l = list(pol)

        if len(c_l) < N: #Si la llista dels seus coeficients es menor que N, s'afegeixen els 0 que calgui al
                              #final.
            zeros = [0]*(N - len(c_l))
            c_l.extend(zeros)

        coef_list.extend(c_l)


    return coef_list


def is_ternary(pol, N, d):

    ternary = True
    index = 0
    num_ones = 0
    num_neg_ones = 0
    num_zeros = 0

    while ternary and index < len(pol):

        if pol[index] == 1:
            num_ones +=1

            if num_ones > d+1:
                ternary = False

        elif pol[index] == -1:
            num_neg_ones +=1

            if num_neg_ones > d:
                ternary = False

        elif pol[index] == 0:
            num_zeros +=1

            if num_zeros > N-d-d:
                ternary = False
        else:
            ternary = False

        index+=1

    return ternary
