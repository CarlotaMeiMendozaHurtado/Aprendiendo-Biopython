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
    
    secuencia_para_blast = registro_descargado.seq

    # --- PASO 4: NCBIWWW (Lanzar un BLAST remoto a los servidores de NCBI) ---
    # Parámetros: programa ("blastn"), base de datos ("nt"), y la secuencia query
    print("Enviando secuencia a BLAST remoto... (puede demorar)")
    resultado_blast_xml = NCBIWWW.qblast("blastn", "nt", secuencia_para_blast)
    
    # --- PASO 5: PARSEAR LOS RESULTADOS DE BLAST ---
    registro_blast = NCBIXML.read(resultado_blast_xml)
    
    # Recorrer los hits encontrados ordenadamente
    for alineamiento in registro_blast.alignments:
        for hsp in alineamiento.hsps:
            # Filtro clásico de examen: Quedarse solo con alineamientos significativos (E-value muy bajo)
            if hsp.evalue < 0.0001:
                print(f"-> Hit encontrado: {alineamiento.title}")
                print(f"   E-value: {hsp.evalue} | Bit-Score: {hsp.bitscore}")