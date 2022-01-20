# Criptosistema NTRU i atacs diversos

Aquest projecte consisteix en la implementació del criptosistema NTRU i els atacs a la vulnerabilitat dels paràmetres, per força bruta, per trobada a mig camí i per reticles al KRP. Està programat amb la versió 9.2 de *SageMath*.

## Fitxers

El projecte disposa d'una versió en format *Script* i d'una altra en format *Notebook*. D'una banda, el *Script* ofereix la possibilitat d'utilitzar el criptosistema i els seus atacs de forma funcional, com qualsevol altre programa. D'altra banda, el *Notebook* pretèn mostrar de forma més didàctica i pedagògica el contingut d'aquest projecte. Amb aquest objectiu, ofereix els resultats d'una execució amb el temps emprat i de comentaris addicionals que en complementen l'explicació.

## Instal·lació i execució

Per tal de poder executar tant el *Script* com el *Notebook* cal instal·lar-se *SageMath* des del link que es mostra a continuació:

```bash
https://www.sagemath.org/

```

# *Notebook*

*Notebook* que implementa el criptosistema NTRU i els atacs a la vulnerabilitat dels paràmetres, per força bruta, per trobada a mig camí i per reticles al KRP. A més a més, ofereix el resultat d'una execució amb el temps emprat i comentaris addicionals que en complementen l'explicació. Està programat amb la versió 9.2 de *SageMath*.

## Execució

Un cop descarregat *SageMath*, des del terminal propi de *SageMath* cal executar la comanda que segueix:

```bash
sage -n jupyter
```

## Fitxers

El *Notebook* es correspon a un únic fitxer on hi apareixen exactament les mateixes classes presents al *Script*, a excepció del codi del fitxer **ntru**.

# *Script*

*Script* que implementa el criptosistema NTRU i els atacs a la vulnerabilitat dels paràmetres, per força bruta, per trobada a mig camí i per reticles al KRP. Està programat amb la versió 9.2 de *SageMath*.

## Execució

Un cop descarregat *SageMath*, des del terminal propi de *SageMath* i dins de la carpeta on es trobi el projecte cal executar la comanda que segueix:

```bash
sage ntru.sage
```

D'aquesta manera es mostraran per pantalla les comandes i els paràmetres necessaris per a cada funcionalitat.

## Fitxers

El *Script* conté un total de 7 fitxers d'extensió *.sage*. Són els següents:

- **ntru**: Conté el main *script*.

- **auxMethods**: Conté diversos mètodes auxiliars que s'empren en les diverses classes del projecte.

- **vulParamAttack**: Fitxer que conté la classe *Param\_vulnerability* referent a l'atac a la vulnerabilitat dels paràmetres.

- **ntruClass**: Fitxer que conté la classe *Ntru* referent al criptosistema NTRU.

- **bfAttack**: Fitxer que conté la classe *Brute\_Force* referent a l'atac per força bruta.

- **mitmAttack**: Fitxer que conté la classe *Meet\_in\_the\_middle* referent a l'atac per trobada a mig camí.

- **krpRetAttack**: Fitxer que conté la classe *Krp\_lattice* referent a l'atac per reticles al KRP.

# Autoria

Aquest projecte ha estat desenvolupat per en Miquel Guiot Cusidó, alumne de la Facultat de Matemàtiques i Informàtica de la Universitat de Barcelona, sota la tutorització del Dr. Xavier Guitart Morales com a part del Treball de Final de Grau.
