load("ntruClass.sage")
load("vulParamAttack.sage")
load("bfAttack.sage")
load("mitmAttack.sage")
load("krpRetAttack.sage")
if __name__ == "__main__":
    text = """NTRU v1.0
        Usage:
          ntru.sage [options] gen_keys N p q d pub_key_file priv_keys_file
          ntru.sage [options] enc pub_key_file message_file enc_message_file
          ntru.sage [options] dec priv_keys_file enc_message_file dec_message_file
          ntru.sage [options] vul_param pub_key_file enc_message_file dec_message_file
          ntru.sage [options] bf pub_key_file enc_message_file dec_message_file
          ntru.sage [options] mitm pub_key_file enc_message_file dec_message_file
          ntru.sage [options] krp_ret pub_key_file enc_message_file dec_message_file
          ntru.sage (-h | --help)
          ntru.sage --version

        Options:
          -h, --help         Show this screen.
          --version          Show version.
        """
    if len(sys.argv) == 1:
        print(text)

    else:

        if len(sys.argv) == 2:

            if sys.argv[1]== "-h" or sys.argv[1]== "--help":
                print(text)

            if sys.argv[1] == "--version":
                print("NTRU 1.0")

        else:

            if sys.argv[1]== "-h" or sys.argv[1]== "--help":
                print(text)

            if sys.argv[1] == "--version":
                print("NTRU 1.0")

            if sys.argv[2]== "-h" or sys.argv[2]== "--help":
                print(text)

            if sys.argv[2] == "--version":
                print("NTRU 1.0")

            if sys.argv[1] == "gen_keys":
                a = Ntru(Integer(sys.argv[2]), Integer(sys.argv[3]), Integer(sys.argv[4]), Integer(sys.argv[5]), True)
                a.gen_keys(sys.argv[6], sys.argv[7])

            if sys.argv[2] == "gen_keys":
                a = Ntru(Integer(sys.argv[3]), Integer(sys.argv[4]), Integer(sys.argv[5]), Integer(sys.argv[6]), True)
                a.gen_keys(sys.argv[7], sys.argv[8])

            if len(sys.argv) > 3:
                if sys.argv[3] == "gen_keys":
                    a = Ntru(Integer(sys.argv[4]), Integer(sys.argv[5]), Integer(sys.argv[6]), Integer(sys.argv[7]), True)
                    a.gen_keys(sys.argv[8], sys.argv[9])

            if sys.argv[1] == "enc":
                a = Ntru(11, 7, 256, 2, False)
                a.set_param_and_pub_key(sys.argv[2])
                a.encrypt_file(sys.argv[3], sys.argv[4])
                print("Fitxer {} encriptat a {}".format(sys.argv[3], sys.argv[4]))

            if sys.argv[2] == "enc":
                a = Ntru(11, 7, 256, 2, False)
                a.set_param_and_pub_key(sys.argv[3])
                a.encrypt_file(sys.argv[4], sys.argv[5])
                print("Fitxer {} encriptat a {}".format(sys.argv[4], sys.argv[5]))

            if len(sys.argv) > 3:
                if sys.argv[3] == "enc":
                    a = Ntru(11, 7, 256, 2, False)
                    a.set_param_and_pub_key(sys.argv[4])
                    a.encrypt_file(sys.argv[5], sys.argv[6])
                    print("Fitxer {} encriptat a {}".format(sys.argv[5], sys.argv[6]))

            if sys.argv[1] == "dec":
                a = Ntru(11, 7, 256, 2, False)
                a.set_param_and_priv_keys(sys.argv[2])
                a.decrypt_file(sys.argv[3], sys.argv[4])
                print("Fitxer {} desencriptat a {}".format(sys.argv[3], sys.argv[4]))

            if sys.argv[2] == "dec":
                a = Ntru(11, 7, 256, 2, False)
                a.set_param_and_priv_keys(sys.argv[3])
                a.decrypt_file(sys.argv[4], sys.argv[5])
                print("Fitxer {} desencriptat a {}".format(sys.argv[4], sys.argv[5]))

            if len(sys.argv) > 3:
                if sys.argv[3] == "dec":
                    a = Ntru(11, 7, 256, 2, False)
                    a.set_param_and_priv_keys(sys.argv[4])
                    a.decrypt_file(sys.argv[5], sys.argv[6])
                    print("Fitxer {} desencriptat a {}".format(sys.argv[5], sys.argv[6]))

            if sys.argv[1] == "vul_param":
                a = Param_vulnerability(sys.argv[2])
                exit = a.attack_file(sys.argv[3], sys.argv[4])
                if exit:
                    print("Fitxer {} atacat per vulnerabilitat de parametres a {}".format(sys.argv[3], sys.argv[4]))

            if sys.argv[2] == "vul_param":
                a = Param_vulnerability(sys.argv[3])
                exit = a.attack_file(sys.argv[4], sys.argv[5])
                if exit:
                    print("Fitxer {} atacat per vulnerabilitat de parametres a {}".format(sys.argv[4], sys.argv[5]))

            if len(sys.argv) > 3:
                if sys.argv[3] == "vul_param":
                    a = Param_vulnerability(sys.argv[4])
                    exit = a.attack_file(sys.argv[5], sys.argv[6])
                    if exit:
                        print("Fitxer {} atacat per vulnerabilitat de parametres a {}".format(sys.argv[5], sys.argv[6]))

            if sys.argv[1] == "bf":
                a = Brute_force(sys.argv[2])
                a.attack_file(sys.argv[3], sys.argv[4])
                print("Fitxer {} atacat per força bruta a {}".format(sys.argv[3], sys.argv[4]))

            if sys.argv[2] == "bf":
                a = Param_vulnerability(sys.argv[3])
                a.attack_file(sys.argv[4], sys.argv[5])
                print("Fitxer {} atacat per força bruta a {}".format(sys.argv[4], sys.argv[5]))

            if len(sys.argv) > 3:
                if sys.argv[3] == "bf":
                    a = Param_vulnerability(sys.argv[4])
                    a.attack_file(sys.argv[5], sys.argv[6])
                    print("Fitxer {} atacat per força bruta a {}".format(sys.argv[5], sys.argv[6]))

            if sys.argv[1] == "mitm":
                a = Meet_in_the_middle(sys.argv[2])
                a.attack_file(sys.argv[3], sys.argv[4])
                print("Fitxer {} atacat per trobada a mig camí a {}".format(sys.argv[3], sys.argv[4]))

            if sys.argv[2] == "mitm":
                a = Param_vulnerability(sys.argv[3])
                a.attack_file(sys.argv[4], sys.argv[5])
                print("Fitxer {} atacat per trobada a mig camí a {}".format(sys.argv[4], sys.argv[5]))

            if len(sys.argv) > 3:
                if sys.argv[3] == "mitm":
                    a = Param_vulnerability(sys.argv[4])
                    a.attack_file(sys.argv[5], sys.argv[6])
                    print("Fitxer {} atacat per trobada a mig camí a {}".format(sys.argv[5], sys.argv[6]))

            if sys.argv[1] == "krp_ret":
                a = Krp_lattice(sys.argv[2])
                a.attack_file(sys.argv[3], sys.argv[4])
                print("Fitxer {} atacat amb reticles al KRP a {}".format(sys.argv[3], sys.argv[4]))

            if sys.argv[2] == "krp_ret":
                a = Krp_lattice(sys.argv[3])
                a.attack_file(sys.argv[4], sys.argv[5])
                print("Fitxer {} atacat amb reticles al KRP a {}".format(sys.argv[4], sys.argv[5]))

            if len(sys.argv) > 3:
                if sys.argv[3] == "krp_ret":
                    a = Krp_lattice(sys.argv[4])
                    a.attack_file(sys.argv[5], sys.argv[6])
                    print("Fitxer {} atacat amb reticles al KRP a {}".format(sys.argv[5], sys.argv[6]))