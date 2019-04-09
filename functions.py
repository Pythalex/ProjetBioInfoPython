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

def update_matrice(ancienne_matrice, utos, plus_proches):
    # n = np.size(ancienne_matrice[0]) - 1
    matrice = np.copy(ancienne_matrice)

    # on retire les lignes et colonnes en trop
    matrice = np.delete(matrice, plus_proches[0], 0)
    matrice = np.delete(matrice, plus_proches[1], 1)
    
    # on met à jour les valeurs des lignes et colonnes fusionnées
    i = plus_proches[0]
    j = plus_proches[1]
    ni = utos[i]
    nj = utos[j]
    for k in range(i+1, np.size(matrice[:, 0])):
        dik = (ancienne_matrice[i, k] if k < j else ancienne_matrice[k, j + 1])
        djk = (ancienne_matrice[k, j] if k < j else ancienne_matrice[k, j + 1])
        matrice[i, k] = (dik * ni + djk * nj) / float(ni + nj)
    for k in range(i - 1, -1, -1):
        dki = ancienne_matrice[k, i]
        dkj = ancienne_matrice[k, j]
        matrice[k, i] = (dki * ni + dkj * nj) / float(ni + nj)
    
    utos[plus_proches[0]] += utos[plus_proches[1]]
    del utos[plus_proches[1]]

    return matrice

def upgma(matrice, dic):
    
    nouvelle = None
    newick = ""
    utos = [1] * np.size(matrice[0])
    while not np.equal(nouvelle, matrice).all():
        i1, i2, d = plus_proche_dans_matrice(matrice)
        print(i1, i2, d)
        newick = "({} : {}, {} : {})".format(dic[i1], d/2, dic[i2], d/2)
        print(newick)
        utos[i1] += utos[i2]
        del utos[i2]
        nouvelle = update_matrice(matrice, utos, (i1, i2))
        print(nouvelle)
        break
    
    return newick

fasta = "cat_dna.fasta"
dic = lire_fasta(fasta)
chiffre_vers_nom = list(dic.keys())
filename = "distmatrix.npy"

if os.path.isfile(filename):
    m = np.load(filename)
else:
    m = matrice_distance(dic)
    save_matrix(filename, m)

#print(m)
#save_matrix(file, m)
# upgma(m, chiffre_vers_nom)
a = np.array([
    [0, 4, 8, 2],
    [0, 0, 8, 4],
    [0, 0, 0, 8],
    [0, 0, 0, 0]
],
dtype=np.float64)
utos = [1, 1, 1, 1, 1]
fusions = []
print("Matrice avant update")
print(a)
print("Matrice après update 1")
fusions.append(plus_proche_dans_matrice(a))
print("plus proches : i {} | j {} | minimum {}".format(*fusions[-1]))
a = update_matrice(a, utos, plus_proche_dans_matrice(a))
print(a)
print("Matrice après update 2")
fusions.append(plus_proche_dans_matrice(a))
print("plus proches : i {} | j {} | minimum {}".format(*fusions[-1]))
a = update_matrice(a, utos, plus_proche_dans_matrice(a))
print(a)
print("Matrice après update 3")
fusions.append(plus_proche_dans_matrice(a))
print("plus proches : i {} | j {} | minimum {}".format(*fusions[-1]))
a = update_matrice(a, utos, plus_proche_dans_matrice(a))
print(a)

print("Fusions : {}".format(fusions))