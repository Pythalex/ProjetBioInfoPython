import re
import numpy as np
import time
import os

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
                    tmp_seq = []

            line = file.readline()

        return dic


class SimilarityMatrix:
    def __init__(self, filename):
        with open(filename) as f:
            self.letters = f.readline().split()
            self.dic = {"A" : 0, "C" : 1, "G" : 2, "T" : 3, "-" : 4}
            self.values = np.loadtxt(filename, skiprows=1, usecols=range(1, len(self.letters) + 1))
        
    def score(self, letter1, letter2): # return the similarity score between letter1 and letter2

        return self.values[self.dic[letter1]][self.dic[letter2]]

D = 0
L = 1
U = 2
O = -1
similarity_matrix = SimilarityMatrix('dna_matrix')

def score_alignement(seq1, seq2):
    global similarity_matrix
    m = similarity_matrix

    if len(seq1) != len(seq2):
        raise Exception("Séquences pas de même longueur")
    
    sum = 0
    for let1, let2 in zip(seq1, seq2):
        sum += m.score(let1, let2)
    return sum


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
    
    if val1 > val2:
        if val1 > val3:
            return val1, 0
        else:
            return val3, 2
    else:
        if val2 > val3:
            return val2, 1
        else:
            return val2, 2

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
    
    return track[::-1]

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

def pdistance(al1, al2):
    nb = 0.
    for e1, e2 in zip(al1, al2):
        nb += 1 if e1 != e2 else 0
    return nb / len(al1)

def juke_cantor(pdistance):
    return -(3/4) * np.log(1 - (4/3) * pdistance)

def matrice_distance(dic):
    keys = list(dic.keys())
    matrix = np.zeros((len(keys), len(keys)), dtype=np.float64)
    print(type(matrix[0,0]))
    
    for i, key1 in enumerate(keys):
        for j in range(i, len(keys)):
            key2 = keys[j]
            al1, al2 = needleman_wunsch(dic[key1], dic[key2])
            pdist = pdistance(al1, al2)
            dist = juke_cantor(pdist)
            matrix[i, j] = float(dist)
            
        print("{}/{}".format(i, len(keys)))
    
    return matrix

def save_matrix(filename, matrix):
    np.save(filename, matrix)

def plus_proche_dans_matrice(matrice):
    min = matrice[0][1]
    
    mini, minj = 0, 1
    for i in range(0, np.size(matrice[:,0])):
        for j in range(i+1, np.size(matrice[0])):
            if min > matrice[i, j]:
                min = matrice[i, j]
                mini = i
                minj = j
    return mini, minj, min

def update_matrice(ancienne_matrice, arbre, plus_proches):
    # n = np.size(ancienne_matrice[0]) - 1
    matrice = np.copy(ancienne_matrice)

    # on retire les lignes et colonnes en trop
    matrice = np.delete(matrice, plus_proches[1], 0)
    matrice = np.delete(matrice, plus_proches[1], 1)
    
    # on met à jour les valeurs des lignes et colonnes fusionnées
    i = plus_proches[0]
    j = plus_proches[1]
    ni = arbre[i].uto()
    nj = arbre[j].uto()
    for k in range(i+1, np.size(matrice[:, 0])):
        dik = (ancienne_matrice[i, k] if k < j else ancienne_matrice[k, j + 1])
        djk = (ancienne_matrice[k, j] if k < j else ancienne_matrice[k, j + 1])
        matrice[i, k] = (dik * ni + djk * nj) / float(ni + nj)
    for k in range(i+1, np.size(matrice[:, 0])):
        dki = ancienne_matrice[i, k]
        dkj = ancienne_matrice[j, k]
        matrice[i, k] = (dki * ni + dkj * nj) / float(ni + nj)

    return matrice

def upgmaOld(matrice, dic):
    
    matrice
    newick = ""
    utos = [1] * np.size(matrice[0])

    lastd = 0

    while not np.size(matrice[0]) == 1:

        i1, i2, d = plus_proche_dans_matrice(matrice)

        newick = "({} : {}, {} : {})".format(dic[i1], d/2 - lastd, dic[i2], d/2)
        lastd = d

        dic[i1] = "({} : {}, {} : {})".format(dic[i1], d/2, dic[i2], d/2)
        del dic[i2]

        matrice = update_matrice(matrice, utos, (i1, i2))

        utos[i1] += utos[i2]
        del utos[i2]
    
    return newick + ";"

class Feuille:
    nom = ""
    distance_gauche = 0
    distance_droite = 0

    def __init__(self, nom):
        self.nom = nom

    def toString(self, first=False):
        return self.nom
    
    def uto(self):
        return 1

class Noeud:
    fils_gauche = None
    fils_droite = None
    distance_gauche = 0
    distance_droite = 0

    def __init__(self, g, dg, d, dd):
        self.fils_gauche = g
        self.fils_droite = d
        self.distance_gauche = dg
        self.distance_droite = dd

    def toString(self, first=True):
        return "({} : {}, {} : {}){}".format(self.fils_gauche.toString(False), 
            self.distance_gauche, 
            self.fils_droite.toString(False), 
            self.distance_droite,
            "" if not first else ";")

    def uto(self):
        return self.fils_gauche.uto() + self.fils_droite.uto()

def upgma(matrice, dic):
    arbre = [Feuille(i) for i in dic]   

    while np.size(matrice[0]) > 1:
        i, j, distance = plus_proche_dans_matrice(matrice)

        matrice = update_matrice(matrice, arbre, (i, j))

        print(arbre[i].toString())
        arbre[i] = Noeud(arbre[i], 
            distance/2 - arbre[i].distance_gauche, 
            arbre[j], 
            distance/2)
        del arbre[j]

        input()


    return arbre[0].toString()

def print_matrix(m):

    print("  ", end="")
    for j in range(np.size(m[0])):
        if j < 10:
            print(" {}  ".format(j), end="")
        else:
            print(" {} ".format(j), end="")
    print()
    for i in range(np.size(m[:, 0])):
        print((" " if i < 10 else "") + str(i), end="")
        print(m[i])

fasta = "cat_dna.fasta"
dic = lire_fasta(fasta)
chiffre_vers_nom = list(dic.keys())
filename = "distmatrix.npy"

if os.path.isfile(filename):
    m = np.load(filename)
else:
    m = matrice_distance(dic)
    save_matrix(filename, m)

# print_matrix(np.around(np.dot(m, 100)))
# n = (upgma(m, chiffre_vers_nom))
# print(n)

# a = np.array([
#     [0, 4, 8, 2],
#     [0, 0, 8, 4],
#     [0, 0, 0, 8],
#     [0, 0, 0, 0]
# ],
# dtype=np.float64)
# utos = [1, 1, 1, 1, 1]
# fusions = []

# print("Matrice avant update")
# print(a)
# print("Matrice après update 1")
# fusions.append(plus_proche_dans_matrice(a))
# print("plus proches : i {} | j {} | minimum {}".format(*fusions[-1]))
# a = update_matrice(a, utos, plus_proche_dans_matrice(a))
# print(a)
# print("Matrice après update 2")
# fusions.append(plus_proche_dans_matrice(a))
# print("plus proches : i {} | j {} | minimum {}".format(*fusions[-1]))
# a = update_matrice(a, utos, plus_proche_dans_matrice(a))
# print(a)
# print("Matrice après update 3")
# fusions.append(plus_proche_dans_matrice(a))
# print("plus proches : i {} | j {} | minimum {}".format(*fusions[-1]))
# a = update_matrice(a, utos, plus_proche_dans_matrice(a))
# print(a)
# print("Fusions : {}".format(fusions))

def get_Dab(matrice, a, b):
    if b < a:
        c = a
        a = b
        b = c
    return matrice[a, b]

def calculer_Sx(matrice, x, arbre):
    return (np.sum(matrice[x]) + np.sum(matrice[:,x])) / (np.size(matrice[0])-2)

def calculer_Mij(matrice, i, j, arbre):
    return get_Dab(matrice, i, j) - calculer_Sx(matrice, i, arbre) - calculer_Sx(matrice, j, arbre)

def plus_petit_Mij(matrice, arbre):
    min = calculer_Mij(matrice, 0, 1, arbre)
    
    mini, minj = 0, 1
    for i in range(0, np.size(matrice[:,0])):
        for j in range(i+1, np.size(matrice[0])):
            mij = calculer_Mij(matrice, i, j, arbre)
            if min > mij:
                min = mij
                mini = i
                minj = j
    return mini, minj, min

def calculer_sau(matrice, a, b, arbre):
    return get_Dab(matrice, a, b) / 2 + (calculer_Sx(matrice, a, arbre) - calculer_Sx(matrice, b, arbre)) / 2

def creer_noeud_ab(matrice, a, b, arbre):
    if b < a:
        c = a
        a = b
        b = c
    
    return Noeud(arbre[a], calculer_sau(matrice, a, b, arbre), arbre[b], calculer_sau(matrice, b, a, arbre))

def cycle(matrice, arbre):
    i, j, min = plus_petit_Mij(matrice, arbre)
    n = creer_noeud_ab(matrice, i, j, arbre)
    arbre[i] = n
    del arbre[j]
    return i, j

def neighbor_joining(matrice, dic):
    arbre = [Feuille(i) for i in dic]

    print(arbre)
    while (np.size(matrice[0])) > 2:
        i, j = cycle(matrice, arbre)
        matrice = update_matrice(matrice, arbre, (i, j))
        print(arbre)
    
    return arbre[0]

def before():
    return (np.array([
        [0, 5, 4, 7, 6, 8],
        [0, 0, 7, 10, 9, 11],
        [0, 0, 0, 7, 6, 8],
        [0, 0, 0, 0, 5, 9],
        [0, 0, 0, 0, 0, 8],
        [0, 0, 0, 0, 0, 0]
    ]
    ),
    [Feuille(i) for i in ["A", "B", "C", "D", "E", "F"]])

def test_calculer_sx():
    m, arbre = before()
    assert calculer_Sx(m, 0, arbre) == 7.5
    assert calculer_Sx(m, 1, arbre) == 10.5
    assert calculer_Sx(m, 2, arbre) == 8
    assert calculer_Sx(m, 3, arbre) == 9.5
    assert calculer_Sx(m, 4, arbre) == 8.5
    assert calculer_Sx(m, 5, arbre) == 11
test_calculer_sx()

def test_calculer_mij():
    m, arbre = before()
    assert calculer_Mij(m, 0, 1, arbre) == -13
    assert calculer_Mij(m, 3, 4, arbre) == -13
test_calculer_mij()

def test_calculer_sau():
    m, arbre = before()
    assert calculer_sau(m, 0, 1, arbre) == 1
    assert calculer_sau(m, 1, 0, arbre) == 4
test_calculer_sau()

def test_creer_noeud_ab():
    m, arbre = before()
    n = creer_noeud_ab(m, 0, 1, arbre)
    assert n.distance_gauche == 1
    assert n.distance_droite == 4
    assert n.fils_gauche.toString() == "A"
    assert n.fils_droite.toString() == "B"
test_creer_noeud_ab()

print(neighbor_joining(m, dic))