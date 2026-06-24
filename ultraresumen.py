''' --------- 1B Operaciones Biológicas con Secuencias ---------'''
from Bio.Seq import Seq

# 1. Crear una secuencia (Ojo: Son inmutables, no puedes hacer mi_adn[0] = 'C')
mi_adn = Seq("ATGCGTACGTGACGTTAA")

# 2. Operaciones biológicas principales
mi_arn = mi_adn.transcribe()            # Reemplaza 'T' por 'U'
adn_vuelta = mi_arn.back_transcribe()   # Reemplaza 'U' por 'T'
comp_inverso = mi_adn.reverse_complement() # Voltea la secuencia y cambia bases (A<->T, C<->G)

# 3. Traducción y sus condiciones especiales (¡Muy preguntable!)
# Por defecto traduce todo e incluye asterisco '*' en el codón de parada (TAA, TAG, TGA)
proteina_normal = mi_adn.translate() 

# Opción A: Detener la traducción al llegar al primer codón de parada (sin incluir el '*')
proteina_limpia = mi_adn.translate(to_stop=True)

# Opción B: Forzar a que sea un CDS estricto (Debe medir múltiplo de 3, empezar por ATG y acabar en parada)
# Si no cumple todas las condiciones biológicas, este parámetro lanzará un error (útil para validar)
proteina_cds = mi_adn.translate(cds=True)



''' --------- 1C Manipulación de Ficheros (SeqIO y SeqRecord) ---------'''
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# 1. Crear un registro estructurado a mano (Secuencia + Metadatos)
registro = SeqRecord(
    Seq("ATGCGTA"), 
    id="ID12345", 
    name="GenEjemplo", 
    description="Secuencia de examen de la UBU"
)

# 2. LEER UN FICHERO: Hay dos funciones obligatorias según el caso
# CASO A: El archivo tiene una ÚNICA secuencia (Si tiene más o ninguna, da error)
registro_unico = SeqIO.read("unasecuencia.fasta", "fasta")
print(registro_unico.id, registro_unico.seq)

# CASO B: El archivo tiene MÚLTIPLES secuencias (Devuelve un iterador/generador)
# Es la opción más segura si no sabes cuántas secuencias vienen
for reg in SeqIO.parse("muchas_secuencias.fasta", "fasta"):
    print(f"ID: {reg.id} | Tamaño: {len(reg.seq)}")

# Si necesitas contar cuántas hay o acceder a una posición concreta, conviértelo a lista:
lista_registros = list(SeqIO.parse("muchas_secuencias.fasta", "fasta"))
primer_registro = lista_registros[0]

# 3. ESCRIBIR/CONVERTIR UN FICHERO
# Toma una lista de registros, el nombre del archivo de salida y el formato
SeqIO.write(lista_registros, "salida.gbk", "genbank")


''' --------- 2B Alineamiento Pareado con Biopython ---------'''
# nw y sw


''' --------- 2D Alineamiento Múltiple ---------'''
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# 1. Crear un objeto MSA a mano (Todas las secuencias deben tener EXACTAMENTE la misma longitud)
alineamiento_multiple = MultipleSeqAlignment([
    SeqRecord(Seq("A-GCT"), id="seq1"),
    SeqRecord(Seq("ATC-T"), id="seq2"),
    SeqRecord(Seq("A-C-T"), id="seq3")
])

# 2. Averiguar cuántas columnas (longitud total con gaps) tiene el alineamiento
columnas = alineamiento_multiple.get_alignment_length()  # Devolverá 5


''' --------- 3B y ENTREZ: BLAST Remoto y Descargas Automáticas de NCBI ---------'''
from Bio import Entrez
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# --- PASO 1: IDENTIFICACIÓN OBLIGATORIA ---
# Si no pones esto, NCBI te bloqueará la IP por seguridad
Entrez.email = "tu_usuario@ubu.es"
Entrez.api_key = "tu_api_key_personal" # Permite subir de 3 a 10 peticiones por segundo

# --- PASO 2: ESEARCH (Buscar identificadores por texto) ---
# Buscaremos en la base de datos de nucleótidos términos relacionados con un organismo
manejador_busqueda = Entrez.esearch(db="nucleotide", term="Cypripedium gansuense", retmax=3)
datos_busqueda = Entrez.read(manejador_busqueda)
manejador_busqueda.close()

lista_ids = datos_busqueda["IdList"] # Contiene una lista de IDs de NCBI (ej. ['123456', '78910'])

# --- PASO 3: EFETCH (Descargar el archivo completo usando los IDs) ---
if lista_ids:
    id_a_descargar = lista_ids[0]
    manejador_descarga = Entrez.efetch(db="nucleotide", id=id_a_descargar, rettype="gb", retmode="text")
    
    # Lo parseamos directamente con SeqIO porque sabemos que viene en formato GenBank ("gb")
    registro_descargado = SeqIO.read(manejador_descarga, "genbank")
    manejador_descarga.close()
    
    queary_secuencia_para_blast = registro_descargado.seq

    # --- PASO 4: NCBIWWW (Lanzar un BLAST remoto a los servidores de NCBI) ---
    # Parámetros: programa ("blastn"), base de datos ("nt"), y la secuencia query
    print("Enviando secuencia a BLAST remoto... (puede demorar)")
    resultado_blast_xml = NCBIWWW.qblast("blastn", "nt", queary_secuencia_para_blast)
    
    # --- PASO 5: PARSEAR LOS RESULTADOS DE BLAST ---
    registro_blast = NCBIXML.read(resultado_blast_xml)
    
    # Recorrer los hits encontrados ordenadamente
    for alineamiento in registro_blast.alignments:
        for hsp in alineamiento.hsps:
            # Filtro clásico de examen: Quedarse solo con alineamientos significativos (E-value muy bajo)
            if hsp.evalue < 0.0001:
                print(f"-> Hit encontrado: {alineamiento.title}")
                print(f"   E-value: {hsp.evalue} | Bit-Score: {hsp.bitscore}")



''' 2 --------- 3B y ENTREZ: BLAST Remoto y Descargas Automáticas de NCBI ---------'''

''' 1. --------- Envia secuencia a internet --------- '''
from Bio.Blast import NCBIWWW

# 1. envía la secuencia a internet
m_m = NCBIWWW.qblast("blastn", "nt" , "ATGCGTACGT")

# 2. lee

datos_xml = m_m.read()
m_m.close

# 3. guardar xml
with open("titulo.xml" , "w") as archivo:
    archivo.write(datos_xml)

''' 2.  --------- Leer Blast, parsear XML --------- '''
from Bio.Blast import NCBIXML
from Bio import Align
from Bio.Align import PairwiseAligner

# 1. Abrir y leer el archivo XML
with open("titulo.xml" , "r") as archivo:
    rec = NCBIXML.read(archivo)

# 2. Agarrar el primer resultado (el mejor Hit)
primer_hit = rec.alingments[0]
print(primer_hit)

# 3. Agarrar la mejor zona alineada (el mejor HSP)
mejor_hsp = primer_hit.hsps[0]
print("E-value:", mejor_hsp.expect)  # Evalúa significancia (lo buscas cercano a 0)
print("Score:", mejor_hsp.score)     # Puntuación de calidad

''' 3.  --------- Descargar de NCBI con Entrez --------- '''
from Bio import Entrez
from Bio import SeqIO

# 1. Identificarse SIEMPRE (Obligatorio para que no te bloqueen)
Entrez.email = "tu_correo@ubu.es"

# 2. Descargar usando el Accession (función efetch)
# Parámetros: base de datos, id (accession), tipo de retorno
m_m = Entrez.efetch(db="nucleotide", id="NC_00123", rettype="gb", retmode="text")

# 3. Convertirlo en un objeto SeqRecord manejable
registro_completo = SeqIO.read(m_m, "genbank")
m_m.close()

# 4. Extraer los datos que necesitas
secuencia_hit = registro_completo.seq  # La secuencia de letras pura
organismo = registro_completo.annotations.get("organism") # El nombre científico

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
from Bio import Entrez
from Bio import SeqIO
from Bio import Align
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# =====================================================================
# PASO 1: CONFIGURAR TU IDENTIFICACIÓN (Obligatorio para NCBI Entrez)
# =====================================================================
Entrez.email = "tu_correo@ubu.es"


# =====================================================================
# PASO 2: EJECUTAR EL BLAST REMOTO Y GUARDAR EL XML
# =====================================================================
mi_secuencia_query = "ATGCGTACGT"  # Tu secuencia problema original

# Enviamos a internet (blastn = nucleótido contra nucleótido, base de datos = nt)
handle_blast = NCBIWWW.qblast("blastn", "nt", mi_secuencia_query)
datos_xml = handle_blast.read()
handle_blast.close()

# Guardamos el XML localmente
with open("resultado_blast.xml", "w") as archivo_xml:
    archivo_xml.write(datos_xml)


# =====================================================================
# PASO 3: PARSEAR EL XML PARA EXTRAER EL MEJOR HIT (ACCESSION)
# =====================================================================
with open("resultado_blast.xml", "r") as archivo_xml:
    blast_record = NCBIXML.read(archivo_xml)

# Agerramos el primer alineamiento (el hit con mayor puntuación)
mejor_hit = blast_record.alignments[0]
accession_mejor_hit = mejor_hit.accession  # Guardamos su identificador único

print(f"El mejor hit encontrado tiene el Accession: {accession_mejor_hit}")


# =====================================================================
# PASO 4: DESCARGAR LA SECUENCIA COMPLETA DESDE NCBI EN_TREZ
# =====================================================================
# Usamos el accession del paso anterior para bajar el archivo GenBank (.gb)
handle_entrez = Entrez.efetch(
    db="nucleotide", 
    id=accession_mejor_hit, 
    rettype="gb", 
    retmode="text"
)

# SeqIO.read lo convierte automáticamente en un objeto estructurado (SeqRecord)
registro_descargado = SeqIO.read(handle_entrez, "genbank")
handle_entrez.close()

secuencia_hit_completa = registro_descargado.seq  # Letras puras del hit
print(f"Organismo descargado: {registro_descargado.annotations.get('organism')}")


# =====================================================================
# PASO 5: ALINEAMIENTO PAREADO GLOBAL (NEEDLEMAN-WUNSCH)
# =====================================================================
# Creamos el alineador automático de Biopython
alineador = Align.PairwiseAligner()

# Elegimos el modo global para que use el algoritmo de Needleman-Wunsch
alineador.mode = "global"  # (Si quisieras Smith-Waterman pondrías 'local')

# Ejecutamos pasando nuestra secuencia original y la que acabamos de descargar
resultados_alineamiento = alineador.align(mi_secuencia_query, secuencia_hit_completa)

# Mostramos el score óptimo final y el dibujo del primer alineamiento
print(f"Puntuación del alineamiento óptimo: {resultados_alineamiento.score}")
print(resultados_alineamiento[0])