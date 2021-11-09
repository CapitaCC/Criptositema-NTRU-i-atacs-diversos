# Criptosistema NTRU

*Script* que implementa el criptosistema NTRU i els atacs a la vulnerabilitat dels paràmetres, per força bruta, per trobada a mig camí i per reticles al KRP. Està programat amb la versió 9.2 de *SageMath*.

## Instal·lació i execució

Per tal de poder executar aquest *script*, en primer lloc cal instal·lar-se *SageMath* des del link que es mostra a continuació:

```bash
https://www.sagemath.org/
```

Acte seguit, des del terminal de *SageMath*i dins de la carpeta on es trobi el projecte cal executar la comanda que segueix:

```bash
sage ntru.sage
```

## Fitxers

El projecte conté un total de 7 fitxers d'extensió .sage. Són els següents:

- **ntru**: Conté el main *script*.

- **auxMethods**: Conté diversos mètodes auxiliars que s'empren en les diverses classes del projecte.

- **vulParamAttack**: Fitxer que conté la classe *Param\_vulnerability* referent a l'atac a la vulnerabilitat dels paràmetres.

- **ntruClass**: Fitxer que conté la classe *Ntru* referent al criptosistema NTRU.

- **bfAttack**: Fitxer que conté la classe *Brute\_Force* referent a l'atac per força bruta.

- **mitmAttack**: Fitxer que conté la classe *Meet\_in\_the\_middle* referent a l'atac per trobada a mig camí.

- **krpRetAttack**: Fitxer que conté la classe *Krp\_lattice* referent a l'atac per reticles al KRP.


## Autoria

Aquest projecte ha estat desenvolupat per en Miquel Guiot Cusidó, alumne de la Facultat de Matemàtiques i Informàtica de la Universitat de Barcelona, sota la tutorització del Dr. Xavier Guitart Morales com a part del Treball de Final de Grau.

