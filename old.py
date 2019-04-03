#%% [markdown]
# # Projet phylogénétique
# Théophile Sanchez (theophile.sanchez@inria.fr) - Sarah Cohen Boulakia
# 
# ------
# 
# ## Introduction
# 
# Au cours de ce projet, vous étudierez trois espèces disparues de félins qui vivaient autrefois sur le continent Américain. Ces trois espèces, le _smilodon_ (tigre à dents de sabre), l'_homotherium_ (_scimitar toothed tigers_) et _M. trumani_ (guépard américain) se sont éteintes il y a environ 13 000 ans, à la fin de la dernière période glaciaire. Des séquences ADN partielles de la protéine cytochrome b de ces espèces ont pu être séquencées et vont vous permettre de retrouver les liens de parentés entre ces espèces et des espèces de félins contemporaines : le chat domestique, le lion, le léopard, le tigre, le puma, le guépard et les chats sauvages africains, chinois et européens. Sont aussi présent dans le jeu de donnée des séquences issues d'espèces extérieures aux félins.
# 
# Afin de reconstruire l'arbre phylogénétique de ces espèces, vous utiliserez une méthode basée sur le calcul des distances évolutives entre les séquences ADN des protéines. Sachez qu'une démarche similaire peut-être appliquée aux séquences d'acides aminés.
# 
# Les différentes étapes qui vous permetterons de construire l'arbre sont détaillées dans ce notebook. Vous devrez implémenter les algorithmes en Python et répondre aux questions dans les cellules prévues.
# 
# Quelques conseils :
# - Utiliser au maximum les fonctions présentes dans les packages de python (sauf si il vous est explicitement demandé de les réimplémenter). Si un problème vous paraît courant, il existe surement déjà une fonction pour le résoudre. Pour ce projet vous serez limité aux packages de base, à Numpy et ETE (seulement pour l'affichage des arbres).
# - Si une partie de votre code ne vous semble pas très explicite, ajoutez des commentaires pour l'expliquer. Une personne qui lit votre code doit pouvoir comprendre son fonctionnement facilement.
# - N'hésitez pas à chercher dans la documentation et sur internet. Cependant, faites attention au plagiat !
# 
# Le projet est à rendre **en binôme** par mail avant le **22/04**. Vous regrouperez votre notebook et les fichiers nécessaires à son fonctionnement dans une archive portant vos noms et prénoms.
#%% [markdown]
# ------
# ## Importation des séquences
# 
# Le format FASTA permet de stocker plusieurs séquences (ADN, ARN ou peptidiques) dans un fichier. Les séquences que vous allez étudier ont été regroupées dans le fichier `cat_dna.fasta`.
# 
# **Exercice 1 :** Écriver une fonction permettant d'importer un fichier au format fasta et de le stocker dans un dictionnaire. Les clés seront les noms des séquences et les valeurs du dictionnaire seront les séquences d'adn.

#%%
import re

import importlib
pb = importlib.util.find_spec("progressbar")
if pb is not None:
    from progressbar import progressbar

def lire_fasta(filename):

    key = re.compile("[>].*")
    sequence = re.compile("[ACTG]+")

    dic = {}
    with open(filename, "r") as file:
        attente_cle = 0
        lecture_sequence = 1
        mode = attente_cle
        line = file.readline()

        current_key = ""
        current_seq = ""

        tmp_seq = []

        while (line != ''):
            if mode == attente_cle:
                if key.match(line):
                    mode = lecture_sequence
                    current_key = line[1:].strip()
            else:
                if sequence.match(line):
                    tmp_seq.append(line.strip())
                else:
                    current_seq = "".join(tmp_seq)
                    mode = attente_cle
                    dic[current_key] = current_seq

            line = file.readline()

        return dic

#%% [markdown]
# ------
# ## Alignement des séquences
# 
# La méthode que vous utiliserez pour calculer l'arbre phylogénétique nécessite de calculer la distance évolutive entre les séquences. Avant de pouvoir les calculer, il faut d'abord aligner les séquences en considérant trois types de mutations :
# - les substitutions (un nucléotide est remplacé par un autre)
# - les insertions
# - les délétions
# Par exemple, les séquences "ACTCCTGA" et "ATCTCGTGA" ont plusieurs alignements possibles : 
# 
# $A_1$ :
# ```
# -ACTCCTGA
# ATCTCGTGA
# ```
# 
# $A_2$ :
# ```
# A-CTCCTGA
# ATCTCGTGA
# ```
# 
# $A_3$ :
# ```
# AC-TCCTGA
# ATCTCGTGA
# ```
# .
# 
# .
# 
# .
# 
# Le "-" désigne un *gap*, c'est à dire un "trou" dans l'alignement qui a été causé par une insertion ou une déletion. On regroupe ces deux types de mutations sous le terme indel.
# 
# Ces alignements correspondent à une multitude d'histoires phylogénétiques différentes. Pour sélectionner le meilleur alignement il faut donc introduire l'hypothèse du maximum de parcimonie qui privilégie l'histoire phylogénétique qui implique le moins d'hypothèses et donc, le moins de changements évolutifs. Par exemple, parmis les trois alignements ci-dessus on preferera l'alignement 2 car il correspond au scénario avec le moins de mutations:
# - l'alignement 1 implique au minimum 1 indel et 3 substitutions
# - l'alignement 2 implique au minimum 1 indel et 2 substitutions
# - l'alignement 3 implique au minimum 1 indel et 3 substitutions
# 
# On peut maintenant définir un score d'identité que l'on va augmenter de 1 lorsque qu'il n'y pas eu de mutation et ainsi obtenir la matrice suivante :
# 
# |   &nbsp;   | A | C | G | T | - |
# |   -   | - | - | - | - | - |
# | **A** | 1 | 0 | 0 | 0 | 0 |
# | **C** | 0 | 1 | 0 | 0 | 0 |
# | **G** | 0 | 0 | 1 | 0 | 0 |
# | **T** | 0 | 0 | 0 | 1 | 0 |
# | **-** | 0 | 0 | 0 | 0 | 0 |
# 
# Cette matrice correspond au modèle d'évolution de l'ADN défini par Jukes et Cantor qui fait l'hypothèse d'un taux de mutation équivalent pour chacun des nucléotides. Cependant, en réalité ces taux ne sont pas les mêmes partout, on sait par exemple que le taux de transition (substitution A$\leftrightarrow$G ou C$\leftrightarrow$T) est différent du taux de transversions (substitution A$\leftrightarrow$T, C$\leftrightarrow$G, C$\leftrightarrow$A ou G$\leftrightarrow$T) et que d'autres facteurs devrait être pris en compte comme la fréquence du nucléotide dans l'ADN. [C'est pour cette raison qu'il existe beaucoup de modèles différents d'écrivant l'évolution de l'ADN.](https://en.wikipedia.org/wiki/Models_of_DNA_evolution) Dans la suite de ce projet nous utiliserons la matrice de similarité $S$ suivante : 
# 
# |   &nbsp;   | A  | C  | G  | T  | -  |
# |   -   | -  | -  | -  | -  | -  |
# | **A** | 10 | -1 | -3 | -4 | -5 |
# | **C** | -1 | 7  | -5 | -3 | -5 |
# | **G** | -3 | -5 | 9  | 0  | -5 |
# | **T** | -4 | -3 | 0  | 8  | -5 |
# | **-** | -5 | -5 | -5 | -5 | -5 |
# 
# **Exercice 2 :** Écriver la fonction permettant de calculer le score entre deux alignements avec la matrice de similarité précédente puis afficher le score des trois alignements $A_1$, $A_2$ et $A_3$. La classe permettant d'importer une matrice et de calculer le score entre deux lettres vous est déjà fournie, la matrice de similarité est stockée dans le fichier `dna_matrix` :
# 

#%%
import numpy as np


class SimilarityMatrix:
    def __init__(self, filename):
        with open(filename) as f:
            self.letters = f.readline().split()
            self.values = np.loadtxt(filename, skiprows=1, usecols=range(1, len(self.letters) + 1))
        
    def score(self, letter1, letter2): # return the similarity score between letter1 and letter2
        return self.values[self.letters.index(letter1)][self.letters.index(letter2)]
    
# Example
similarity_matrix = SimilarityMatrix('dna_matrix')

#%%
def score_alignement(seq1, seq2):
    
    if len(seq1) != len(seq2):
        raise Exception("Séquences pas de même longueur")
    
    m = SimilarityMatrix('dna_matrix')
    sum = 0
    for let1, let2 in zip(seq1, seq2):
        sum += m.score(let1, let2)
    return sum


#%% [markdown]
# ------
# ### Algorithme de Needleman-Wunsch
# 
# Maintenant que vous avez vu ce qu'est une matrice de similarité et comment calculer le score de similarité d'un alignement, vous allez devoir implémenter un algorithme permettant de trouver le meilleur alignement global entre deux séquences. Avec deux séquences à aligner de taille $n$ et $m$, la première étape consiste à initialiser deux matrices de taille $n \times m$. La première est la matrice de score $M$ et la seconde sera la matrice de *traceback* $T$. 
# 
# Par exemple, avec la matrice $S$ et les séquences $A =$ "ACTCCTGA" et $B =$ "ATCTCGTGA", on initialise $M$ comme si l'on ajoutait des *gaps* partout :
# 
# |   &nbsp;   | - | A | T | C | T | C | G | T | G | A |
# |   -   | - | - | - | - | - | - | - | - | - | - |
# | **-** | 0 |-5 |-10|-15|-20|-25|-30|-35|-40|-45|
# | **A** |-5 | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; |
# | **C** |-10| &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; |
# | **T** |-15| &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; |
# | **C** |-20| &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; |
# | **C** |-25| &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; |
# | **T** |-30| &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; |
# | **G** |-45| &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; |
# | **A** |-40| &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; ||
# 
# Puis on initialise $T$ :
# 
# |   &nbsp;   | - | A | T | C | T | C | G | T | G | A |
# |   -   | - | - | - | - | - | - | - | - | - | - |
# | **-** | o | l | l | l | l | l | l | l | l | l |
# | **A** | u | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; |
# | **C** | u | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; |
# | **T** | u | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; |
# | **C** | u | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; |
# | **C** | u | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; |
# | **T** | u | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; |
# | **G** | u | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; |
# | **A** | u | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; ||
# 
# 
# Il faut ensuite remplir la matrice $M$ en suivant la formule $M_{ij} = \max(M_{i-1j-1} + s(A_i, B_j), M_{ij-1} + s(A_i, gap), M_{i-1j} + s(B_j,gap) )$ avec $i \in {2, \dots, n}$, $j \in {2, \dots, m}$ et la fonction $s$ qui calcule le score de similarité entre deux nucléotides. Pour chaque case de $T$ on remplie par :
# - 'd' (*diagonal*) si $M_{ij}$ a été calculé en utilisant la diagonale $M_{i-1j-1}$,
# - 'l' (*left*) si $M_{ij}$ a été calculé en utilisant la case de gauche $M_{ij-1}$,
# - 'u' (*up*) si $M_{ij}$ a été calculé en utilisant la case du haut $M_{i-1j}$.
# 
# On obtient alors les matrices suivantes $M$ et $T$ : 
# 
# |   &nbsp;   | - | A | T | C | T | C | G | T | G | A |
# |   -   | - | - | - | - | - | - | - | - | - | - |
# | **-** |  0| -5|-10|-15|-20|-25|-30|-35|-40|-45|
# | **A** | -5| 10|  5|  0| -5|-10|-15|-20|-25|-30|
# | **C** |-10|  5|  7| 12|  7|  2| -3| -8|-13|-18|
# | **T** |-15|  0| 13|  8| 20| 15| 10|  5|  0| -5|
# | **C** |-20| -5|  8| 20| 15| 27| 22| 17| 12|  7|
# | **C** |-25|-10|  3| 15| 17| 22| 22| 19| 14| 11|
# | **T** |-30|-15| -2| 10| 23| 18| 22| 30| 25| 20|
# | **G** |-35|-20| -7|  5| 18| 18| 27| 25| 39| 34|
# | **A** |-40|-25|-12|  0| 13| 17| 22| 23| 34| 49|
# 
# |   &nbsp;   | - | A | T | C | T | C | G | T | G | A |
# |   -   | - | - | - | - | - | - | - | - | - | - |
# | **-** | o | l | l | l | l | l | l | l | l | l |
# | **A** | u | d | l | l | l | l | l | l | l | d |
# | **C** | u | u | d | d | l | d | l | l | l | l |
# | **T** | u | u | d | l | d | l | l | d | l | l |
# | **C** | u | u | u | d | l | d | l | l | l | l |
# | **C** | u | u | u | d | d | d | d | d | l | d |
# | **T** | u | u | d | u | d | l | d | d | l | l |
# | **G** | u | u | u | u | u | d | d | u | d | l |
# | **A** | u | d | u | u | u | d | u | d | u | d |
# 
# Il suffit maintenant de regarder le dernier élément $M_{nm} = 49$ pour avoir le score de l'alignement. Pour avoir l'alignement lui-même, il faut partir de $T_{nm}$ et remonter la "trace" jusqu'à arriver au 'o'. Un 'd' correspond à un *match* entre les deux séquences, 'l' à un *gap* dans la séquence $A$ et 'u' à un *gap* dans la séquence $B$. En revenant à l'exemple précédent on obtient la trace suivante :
# 
# |   &nbsp;   | - | A | T | C | T | C | G | T | G | A |
# |   -   | - | - | - | - | - | - | - | - | - | - |
# | **-** | o | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; |
# | **A** | &nbsp; | d | l | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; |
# | **C** | &nbsp; | &nbsp; | &nbsp; | d | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; |
# | **T** | &nbsp; | &nbsp; | &nbsp; | &nbsp; | d | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; |
# | **C** | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | d | &nbsp; | &nbsp; | &nbsp; | &nbsp; |
# | **C** | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | d | &nbsp; | &nbsp; | &nbsp; |
# | **T** | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | d | &nbsp; | &nbsp; |
# | **G** | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | d | &nbsp; |
# | **A** | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | &nbsp; | d |
# 
# Elle correspond à l'alignement :
# ```
# A-CTCCTGA
# ATCTCGTGA
# ```
# 
# **Exercice 3 :** Implémenter l'algorithme de Needlman et Wunsch. Il prendra en paramètre deux séquences et une matrice de similarité et retournera leur alignement. Tester le avec les séquences "ACTCCTGA" et "ATCTCGTGA".

#%%
D = 0
L = 1
U = 2
O = -1

def init_gap(matrix):
    matrix[0] = np.array([-i * 5 for i in range(0, np.size(matrix[0]))])
    matrix[:,0] = np.array([-i * 5 for i in range(0, np.size(matrix[:,0]))])

def init_t(T):
    global L, D, U, O
    
    T[0] = np.array([L]*np.size(T[0]))
    T[:,0] = np.array([U]*np.size(T[:,0]))
    T[0,0] = O
    
def formule_needleman(i, j, M, seq1, seq2):
    global similarity_matrix
    global L, D, U
    s = similarity_matrix.score
    
    val1 = M[i-1, j-1] + s(seq1[i - 1], seq2[j - 1])
    
    val2 = M[i, j-1] + s(seq1[i - 1], '-')
    
    val3 = M[i-1, j] + s(seq2[j - 1], '-')
    
    values = [val1, val2, val3]
    argmax = np.argmax(values)
    return values[argmax], argmax

def backtrack(T):
    global L, D, U, O
    
    i, j = np.size(T[:, 0]) - 1, np.size(T[0]) - 1
    track = []
    
    while T[i, j] != O:
        track.append(T[i, j])
        if T[i, j] == D:
            i, j = i-1, j-1
        elif T[i, j] == L:
            j -= 1
        elif T[i, j] == U:
            i -= 1
        else:
            print("bug")
    
    return np.flip(track)

def alignement_from_backtrack(seq1, seq2, backtrack):
    global D, L, U
    res1, res2 = [], []
    
    c1, c2 = 0, 0
    for i, deplacement in enumerate(backtrack):
        if deplacement == D:
            res1.append(seq1[c1])
            res2.append(seq2[c2])
            c1 += 1
            c2 += 1
        elif deplacement == L:
            res1.append('-')
            res2.append(seq2[c2])
            c2 += 1
        elif deplacement == U:
            res1.append(seq1[c1])
            res2.append('-')
            c1 += 1
        else:
            print("bug")
            
    return "".join(res1), "".join(res2)
            

def needleman_wunsch(seq1, seq2):
    global L, D, U, O
    
    M = np.zeros(shape=(len(seq1) + 1, len(seq2) + 1))
    init_gap(M)
    
    T = np.full((len(seq1) + 1, len(seq2) + 1), O)
    init_t(T)
    
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            M[i, j], T[i, j] = formule_needleman(i, j, M, seq1, seq2)
    
    backt = backtrack(T)
    
    alignement = alignement_from_backtrack(seq1, seq2, backt)
    
    return alignement

#%% [markdown]
# ----
# ## Matrice de distance
# 
# Dans le cas de séquences très proches, on estime que la distance évolutive réelle entre les séquences est proche de la p-distance qui est simplement le nombre de substitution dans l'alignement sur le nombre total de nucléotide. Pour simplifier, on ignore les positions alignées à des gaps. On applique ensuite la correction de Jukes-Cantor afin de prendre en compte le phénomène de saturation (un même site peut muter plusieurs fois au cours du temps). Sa formule est $-(\frac{3}{4})\ln(1-(\frac{4}{3})\times \textit{p-distance})$.
# 
# **Exercice 4 :** Implémenter la fonction retournant la matrice de distance à partir d'un dictionnaire de séquences. 

#%%
def pdistance(al1, al2):
    nb = 0.
    for e1, e2 in zip(al1, al2):
        nb += 1 if e1 != e2 else 0
    return nb / len(al1)

def juke_cantor(pdistance):
    return -(3/4) * np.log(1 - (4/3) * pdistance)

def matrice_distance(dic):
    
    keys = list(dic.keys())

    matrix = np.zeros((len(keys), len(keys)))
    
    if not pb:
        for i, key1 in enumerate(keys):
            for j, key2 in enumerate(keys[i+1:]):
                al1, al2 = needleman_wunsch(dic[key1], dic[key2])
                pdist = pdistance(al1, al2)
                dist = juke_cantor(pdist)
                matrix[i, i+j] = dist
    else:
        for i in progressbar(range(len(keys))):
            key1 = keys[i]
            for j, key2 in enumerate(keys[i+1:]):
                al1, al2 = needleman_wunsch(dic[key1], dic[key2])
                pdist = pdistance(al1, al2)
                dist = juke_cantor(pdist)
                matrix[i, i+j] = dist
            
                
    
    return matrix
    
if __name__ == '__main__':

    print(matrice_distance(lire_fasta("petit.fasta")))

#%% [markdown]
# ------
# ## Construction d'un arbre avec UPGMA
# 
# Grâce aux mesures de distances entre les séquences, on peut maintenant de construire l'arbre phylogénétique des globines. Vous allez devoir pour cela implémenter l'algorithme UPGMA (*unweighted pair group method with arithmetic mean*) qui, malgré son nom compliqué, est l'une des méthodes les plus simples pour la construction d'arbre.
# 
# ### Le format Newick
# 
# Le format Newick est l'un des formats utilisé en phylogénie pour représenter un arbre sous la forme d'une chaine de caractère. Le principe est simple, les groupes ayant la même racine sont écrit entre parenthèses et séparés par des virgules. Un groupe peut être soit une feuille de l'arbre (dans notre cas une séquence), soit un autre groupe. La longueur de la branche de chaque groupe est écrite après un double point et l'arbre est terminé par un point virgule. Pour afficher l'arbre on peut utiliser les fonction du package ETE : 

# #%%
# from ete3 import Tree, TreeStyle

# newick_tree = '((A:2,(B:2,C:3):5):5,D:4);'
# t = Tree(newick_tree)
# ts = TreeStyle()
# ts.show_branch_length = True
# t.render('%%inline', w=183, units='mm', tree_style=ts)

#%% [markdown]
# **Exercice 5 :** Reécriver l'arbre suivant au format Newick puis afficher-le. Les nombres correspondent aux longueurs des branches :
# ![](tree.png)

#%%
# À remplir

#%% [markdown]
# **Exercice 6 :** Expliquer la relation de parenté entre $A$, $B$ et $C$ ? Qu'est ce qui pourrait expliquer ce type d'embranchement dans un arbre ? Donner une réponse détaillée.
#%% [markdown]
# Réponse : 
#%% [markdown]
# ### UPGMA
# 
# L'algorithme UPGMA se base sur la matrice de distance entre les séquences. À chaque itération, les séquences avec la distance la plus faible sont regroupées puis une nouvelle matrice de distance est calculée avec le nouveau groupe. Cette étape est répétée jusqu'à n'avoir plus qu'un seul groupe. Par exemple, avec la matrice de distance entre les séquences $A$, $B$, $C$ et $D$ suivante :
# 
# |   &nbsp;   | A | B | C | D |
# |   -   | - | - | - | - |
# | **A** | &nbsp; | &nbsp; | &nbsp; | &nbsp; |
# | **B** | 4 | &nbsp; | &nbsp; | &nbsp; |
# | **C** | 8 | 8 | &nbsp; | &nbsp; |
# | **D** | 2 | 4 | 8 | &nbsp; |
# 
# Les séquences $A$ et $D$ sont les plus proches ($distance(A,D)=1$). On les regroupe et on met à jour la matrice :
# 
# |   &nbsp;   | (A, D) | B | C |
# |   -   | - | - | - |
# | **(A, D)** | &nbsp; | &nbsp; | &nbsp; |
# | **B** | 4 | &nbsp; | &nbsp; |
# | **C** | 8 | 8 | &nbsp; | &nbsp; |
# 
# On regroupe maintenant $(A,D)$ et $B$ (distance((A,D),B) = 2) :
# 
# |   &nbsp;   | ((A, D), B) | C |
# |   -   | - | - |
# | **((A, D), B)** | &nbsp; | &nbsp; |
# | **C** | 8 | &nbsp; |
# 
# Important : les nouvelles distances sont calculées en moyennant les distances entre les membres du nouveau groupe et des groupes non modifiés pondéré par le nombre d'UTOs dans chaque groupe. Avec $i$ et $j$ les deux membres du groupe nouvellement formé et k les groupes restant : $d_{ij,k} = \frac{n_id_{ik}}{n_i + n_j}+ \frac{n_jd_{jk}}{n_i + n_j}$. Par exemple avec la distance entre $((A, D), B)$ et $C$:
# 
# $distance(((A, D), B), C) = (distance((A, D), C)*2 + distance(B, C)) \mathbin{/} 3 = (8*2 + 8) \mathbin{/} 3 = 8 $.
# 
# L'arbre final écrit dans le format Newick est : $((A, D), B), C);$ 
# 
# Et avec les distances : $((A:1, D:1):1, B:2):2, C:4);$ 
# 
# **Exercice 7 :** Implémenter une version d'UPGMA qui calcule l'arbre au format Newick **avec les distances** puis appliquer votre algorithme aux données. 

#%%
# À remplir

#%% [markdown]
# **Exercice 8 :** Quelles sont les hypothèses faites par UPGMA ? Semblent-elles respectées dans le cas présent ?
#%% [markdown]
# Réponse : 
#%% [markdown]
# ----
# ## Enracinement de l'arbre
# 
# Après avoir utilisé UPGMA pour réaliser votre arbre, l'enracinement s'est normalement fait au poids moyen. 
# 
# **Exercice 9 :** Quelle autre méthode est-il possible d'utiliser pour enraciner un arbre ? Pouvez-vous l'utilisez ici ? Si oui, afficher le nouvel arbre.
#%% [markdown]
# Réponse : 

#%%
# À remplir

#%% [markdown]
# ----
# ## Neighbor-joining
# 
# Le neighbor-joining est un autre algorithme permettant de calculer un arbre phylogénique à partir d'une matrice de distance. Il a l'avantage de faire moins d'hypothèse qu'UPGMA sur les données (elles ne sont plus forcément ultramétrique) et il donne donc de meilleurs arbres dans presque tout les cas. Vous trouverez un example d'application de cet algorithme [ici](http://www.evolution-textbook.org/content/free/tables/Ch_27/T11_EVOW_Ch27.pdf).
# 
# **Exercice 10 :** Implémenter l'algorithme du neighbor-joining et appliquer-le aux données.

#%%
# À remplir

#%% [markdown]
# ----
# ## Conclusion
# 
# **Exercice 11 :** Quelles sont vos conclusions par rapport à l'arbre phylogénique de _smilodon_, _homotherium_ et _M. trumani_ ? Comparer les deux méthodes. Comment expliquer les caractéristiques morphologiques similaires entre ces trois espèces ? Une réponse détaillée est attendue.
#%% [markdown]
# Réponse :

#%%



