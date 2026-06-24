from Bio.Seq import Seq
from Bio import Align
from Bio.Align import substitution_matrices

import numpy as np

''' 
arr = np.array([1,2,3])
sol = arr + 2
#[3,4,5]
''' 




def smith_waterman(self, seq1, seq2, matrix_name, gap_penalty):
    ''' 
    ----------- Preparación -----------
    1. Matriz de sustitucion: que usaremos en este caso, se supone que SW es local, por lo tanto habría que usar Blosum
    2. Para poder crear la matriz, la seq debe ser un str y debemos saber su tamaño para el tamaño del tablero
    3. Creamos las dos matrices llenas de ceros
    4. variables vacias a las que iremos metiendo segun recorramos las matrices
    '''
    # 1. -----------
    matrix = substitution_matrices.load(matrix_name)
    
    # 2. -----------
    s1 = str(seq1)
    s2 = str(seq2)
    # s1 , s2 = str(seq1) , str(seq2)
    t1 = len(s1)
    t2 = len(s2)
    #t1 , t2 = len(t1), len(t2)
    # Basicamente t1 y t2 son las filas y columnas

    # 3. -----------
    m_puntos = np.zeros((t1+1),(t2+1))
    m_pointer = np.zeros((t1+1),(t2+1))
    # matrices las llenamos de ceros y las añadimos 1 f y c mas
    # matrices m_puntos y m_pointer
    # casilla m_puntos[i][j]

    # 4. -----------
    max_score = 0 
    # variable para ir metiendo puntuacion
    max_pos = (0,0)
    # el final del recorrido y donde está el puntaje max

    '''
    ----------- Rellenar puntos ----------
    1. range( 1 , t1+1 ): el 1 evita empezar en el - inicial
    2. match -> diagonal | delete -> cielo (x suprimir, gap en s2) | insert -> isquierda (vertical intro gap en s1)

    1. recorrer las filas y columnas sin tocar lo de ceros y recorriendo eso
    2. vemos sus puntuaciones en match, delete o insert, eligiendo la meno diferente a -x
    3. SI LAS TRES NEGATIVAS, PONER UN 0
    4. Poner el valor
    '''
    '''
    ----------- Rellenar pointer ----------
    5. Condicionales de direccion

        Si 0        -> guarda   0   stop
        Si match    -> guarda   1   diagonal
        Si delete   -> guarda   2   arrriba
        Si insert   -> guarda   3   izquierda
    '''

    # 1.
    for i in range(1 , t1+1):
        for j in range(1,t2+1):
            '''
            match -> significa `nos importan valores que vengan en diagonal`
            m_puntos[i-1][j-1] -> hay algun valor en la daigonal anterior? 
            matrix[s1[i-1]][s2[j-1]] -> sacar el valor que es tomar esa diagonal acorde a la matriz
            (esto está puesto porque a papel lo entiendo pero a codigo no)
            '''
            # 2.
            match = m_puntos[i-1][j-1] + matrix[s1[i-1]][s2[j-1]]
            delete = m_puntos[i-1][j] + gap_penalty
            insert = m_puntos[i][j-1] + gap_penalty

            # 3. coincidion ORO,devolvera el mayor
            mayor = max(0,match, delete, insert)

            # 4. asignar valor a esa posicion en la que estamos
            m_puntos[i][j] = mayor


            # 5. DIRECCION
            if mayor == 0:
                m_pointer[i][j] = 0
            elif mayor == match:
                m_pointer[i][j] = 1
            elif mayor == delete:
                m_pointer[i][j] = 2
            elif mayor == insert:
                m_pointer[i][j] = 3

            # 6. Actualizar variables
            if mayor > max_score:
                max_score = mayor
                max_pos = (i,j)
            # se registra y se memoriza posicion 
    
    '''
    Backtracking
    1. Dos str vacios
    '''

    ali1 , ali2 = "" , "" # guarda ali final

    # empezamos obligatoriamente en max posicion obtenida
    i , j = max_pos # transporta maximo

    # si flecha es 0, paramos

    while m_pointer[i][j] !=0:

        # diagonal
        if m_pointer[i][j] == 1:
            ali1 = s1[i-1] + ali1
            ali2 = s2[j-1] + ali2
            i -=1
            j -=1

        # arriba
        if m_pointer[i][j] == 2:
            ali1 = s1[i-1] + ali1
            ali2 = "-" + ali2
            i -=1

        # izquierda
        if m_pointer[i][j] == 1:
            ali1 = "-" + ali1
            ali2 = s2[j-1] + ali2
            j -=1

        return ali1 , ali2, max_score

            





from Bio import Align
from Bio.Align import substitution_matrices


def needleman_wunch(self, seq1, seq2, matrix_name, gap_penalty):
    # 1.
    matrix = substitution_matrices.load(matrix_name)

    # 2.
    sa , sb = str(seq1) , str(seq2)
    la , lb = len(sa) , len(sb)         # tamaños

    # 3.
    puntuacion = np.zeros((la+1),(lb+1))
    flechas = np.zeros((la+1),(lb+1))

    # 4. no necesarias

    # Recorrer inicial 
    # fila
    for i in range(1,la+1):
        puntuacion[i][0] = puntuacion[i-1][0] + gap_penalty
        flechas[i][0] = 2 # porque estás moviendote hacia abajo

    # columna
    for j in range (1,lb+1):
        puntuacion[0][j] = puntuacion[0][j-1] + gap_penalty
        flechas[0][j] = 3 # porque vienen de la izuiqerda

    # general
    for i in range(1,la+1):
        for j in range(1,lb+1):
            match = puntuacion[i-1][j-1] + matrix[sa[i-1]][sb[j-1]]
            delete = puntuacion[i-1][j-1] + gap_penalty
            insert = puntuacion[i-1][j-1] + gap_penalty

            mejor = max(match,delete,insert)
            puntuacion[i][j] = mejor

    
            if mejor == match:
                flechas[i][j] = 1   
            elif mejor == delete:
                flechas[i][j] = 2   
            else:
                flechas[i][j] = 3  

    # --- BACKTRACKING NEEDLEMAN-WUNSCH ---
    align1, align2 = "", ""

    # Regla 1: Empezamos estrictamente en la esquina inferior derecha
    i, j = la, lb  

    # Regla 2: No paramos hasta que AMBOS índices lleguen a la esquina superior izquierda (0,0)
    while i > 0 or j > 0:
        
        # Caso ideal: Aún nos quedan letras en ambas secuencias para evaluar la flecha
        if i > 0 and j > 0:
            if flechas[i][j] == 1:     # Diagonal
                align1 = sa[i-1] + align1
                align2 = sb[j-1] + align2
                i -= 1
                j -= 1

            elif flechas[i][j] == 2:   # Arriba
                align1 = sa[i-1] + align1
                align2 = "-" + align2
                i -= 1

            elif flechas[i][j] == 3:   # Izquierda
                align1 = "-" + align1
                align2 = sb[j-1] + align2
                j -= 1
                
        # Caso extremo 1: j llegó a 0 (se acabaron las letras de sb), solo podemos subir metiendo gaps en align2
        elif i > 0:
            align1 = sa[i-1] + align1
            align2 = "-" + align2
            i -= 1
            
        # Caso extremo 2: i llegó a 0 (se acabaron las letras de sa), solo podemos ir a la izquierda metiendo gaps en align1
        elif j > 0:
            align1 = "-" + align1
            align2 = sb[j-1] + align2
            j -= 1

    # La puntuación final de NW está guardada en la última celda de la matriz
    puntuacion_final = puntuacion[la][lb]

    return align1, align2, puntuacion_final





# ----------------------------------------------------------------------------------------------------------------------
def smith_waterman(self, seq1, seq2, matrix_name, gap_penalty):
    '''----------- Inicializacion ----------'''
    # 1. -----------
    matrix = substitution_matrices.load(matrix_name)
    # 2. -----------
    s1 , s2 = str(seq1) , str(seq2)
    t1 , t2 = len(t1), len(t2)
    # 3. -----------
    m_puntos = np.zeros((t1+1),(t2+1))
    m_pointer = np.zeros((t1+1),(t2+1))
    # 4. -----------
    max_score = 0 
    # variable para ir metiendo puntuacion
    max_pos = (0,0)
    # el final del recorrido y donde está el puntaje max

    '''----------- Rellenar puntos ----------'''
    # 1.
    for i in range(1 , t1+1):
        for j in range(1,t2+1):
            # 2.
            match = m_puntos[i-1][j-1] + matrix[s1[i-1]][s2[j-1]]
            delete = m_puntos[i-1][j] + gap_penalty
            insert = m_puntos[i][j-1] + gap_penalty
            # 3. coincidion ORO,devolvera el mayor
            mayor = max(0,match, delete, insert)
            # 4. asignar valor a esa posicion en la que estamos
            m_puntos[i][j] = mayor
            # 5. DIRECCION
            if mayor == 0:
                m_pointer[i][j] = 0
            elif mayor == match:
                m_pointer[i][j] = 1
            elif mayor == delete:
                m_pointer[i][j] = 2
            elif mayor == insert:
                m_pointer[i][j] = 3
            # 6. Actualizar variables # se registra y se memoriza posicion 
            if mayor > max_score:
                max_score = mayor
                max_pos = (i,j) 
    
    '''----------- Backtracking ----------'''
    ali1 , ali2 = "" , "" # guarda ali final
    i , j = max_pos # transporta maximo # empezamos obligatoriamente en max posicion obtenida

    # si flecha es 0, paramos
    while m_pointer[i][j] !=0:
        # diagonal
        if m_pointer[i][j] == 1:
            ali1 = s1[i-1] + ali1
            ali2 = s2[j-1] + ali2
            i -=1
            j -=1
        # arriba
        if m_pointer[i][j] == 2:
            ali1 = s1[i-1] + ali1
            ali2 = "-" + ali2
            i -=1
        # izquierda
        if m_pointer[i][j] == 1:
            ali1 = "-" + ali1
            ali2 = s2[j-1] + ali2
            j -=1
        return ali1 , ali2, max_score



# ----------------------------------------------------------------------------------------------------------------------
def needleman_wunch(self, seq1, seq2, matrix_name, gap_penalty):
    ''' Inicializar '''
    # 1.
    matrix = substitution_matrices.load(matrix_name)
    # 2.
    sa , sb = str(seq1) , str(seq2)
    la , lb = len(sa) , len(sb)         # tamaños
    # 3.
    puntuacion = np.zeros((la+1),(lb+1))
    flechas = np.zeros((la+1),(lb+1))
    # 4. no necesarias

    ''' Rellenar '''
    # fila
    for i in range(1,la+1):
        puntuacion[i][0] = puntuacion[i-1][0] + gap_penalty
        flechas[i][0] = 2 # porque estás moviendote hacia abajo

    # columna
    for j in range (1,lb+1):
        puntuacion[0][j] = puntuacion[0][j-1] + gap_penalty
        flechas[0][j] = 3 # porque vienen de la izuiqerda

    # general
    for i in range(1,la+1):
        for j in range(1,lb+1):
            match = puntuacion[i-1][j-1] + matrix[sa[i-1]][sb[j-1]]
            delete = puntuacion[i-1][j-1] + gap_penalty
            insert = puntuacion[i-1][j-1] + gap_penalty
            mejor = max(match,delete,insert)
            puntuacion[i][j] = mejor
            if mejor == match:
                flechas[i][j] = 1   
            elif mejor == delete:
                flechas[i][j] = 2   
            else:
                flechas[i][j] = 3  

    ''' Bcakteacking '''
    align1, align2 = "", ""

    # Regla 1: Empezamos estrictamente en la esquina inferior derecha
    i, j = la, lb  
    # Regla 2: No paramos hasta que AMBOS índices lleguen a la esquina superior izquierda (0,0)
    while i > 0 or j > 0:
        
        # Caso ideal: Aún nos quedan letras en ambas secuencias para evaluar la flecha
        if i > 0 and j > 0:
            if flechas[i][j] == 1:     # Diagonal
                align1 = sa[i-1] + align1
                align2 = sb[j-1] + align2
                i -= 1
                j -= 1
            elif flechas[i][j] == 2:   # Arriba
                align1 = sa[i-1] + align1
                align2 = "-" + align2
                i -= 1
            elif flechas[i][j] == 3:   # Izquierda
                align1 = "-" + align1
                align2 = sb[j-1] + align2
                j -= 1
        # Caso extremo 1: j llegó a 0 (se acabaron las letras de sb), solo podemos subir metiendo gaps en align2
        elif i > 0:
            align1 = sa[i-1] + align1
            align2 = "-" + align2
            i -= 1
        # Caso extremo 2: i llegó a 0 (se acabaron las letras de sa), solo podemos ir a la izquierda metiendo gaps en align1
        elif j > 0:
            align1 = "-" + align1
            align2 = sb[j-1] + align2
            j -= 1

    # La puntuación final de NW está guardada en la última celda de la matriz
    puntuacion_final = puntuacion[la][lb]
    return align1, align2, puntuacion_final
