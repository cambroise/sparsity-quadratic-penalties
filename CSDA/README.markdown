
# Quadrupen


## Introduction

- Interprétation géométrique et aspect algorithmique, pourquoi c’est bien:
    - c’est générique
    - c’est  précis
    - dans quel cadre c’est bien
        - problèmes de taille intermédiaire

## Arguments

- résolution exacte
    - intérêt pour l’interprétabilité
    - voir la conclusion
- point de vue unifiant
- dualité
    - qui permet le calcul de la borne
    - on remplace le sous gradient par le worst case
- créativité (attention surtout valable en 2d
- stabilité
- linéaire par morceau
- ...
- code
- expériences
- borne
    - mais moins bonne que Fenchel
- plus de sous gradient

## Contre-Arguments

- rien de vraiment neuf
- coût quadratique

## Conclusion

- c’est un package qui devrait être utilisé pour interpréter les variables sélectionnées
    - parcque c’est plus précis
    - parce qu’en conséquences on gagne en précision sur le support, mais pas en terme de MSE

## Consigne CSDA

- The journal requires that the manuscript contain either computational or data analysis component.  Papers which are purely theoretical are not appropriate.  Manuscripts describing simulation studies must
    - 1. be thorough with regard to the choice of parameter settings,
    - 2. not over-generalize the conclusions,
    - 3. carefully describe the limitations of the simulation studies, and
    - 4. should guide the user regarding when the recommended methods areappropriated.It is recommended that the author indicate why comparisons cannot be made theoretically.

## TODO

- Reprendre l’intro
    - une vue unifiante
    - fondée sur une dualité alternative
    - qui permet
        - une vision géométrique
            - créativité
        - une analyse sans sous gradient
        - solution exacte
            - précision
            - stabilité
            - le calcul d’une borne
    - analyse de la vitesse de convergence difficile
        - donc nous faisons des simulations
