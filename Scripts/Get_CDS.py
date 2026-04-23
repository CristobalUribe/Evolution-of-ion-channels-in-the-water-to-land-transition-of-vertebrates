import os
import argparse
from Bio import SeqIO
import re


# Dictionary of species names to CDS files
MAPA_ARCHIVOS = {
    "Aguila": "Aguila_CDS", "Cisnevul": "Cisnevul_CDS", "Gallo": "Gallo_CDS", "Mardo": "Mardo_CDS", "Pavo": "Pavo_CDS", "Takfla": "Takfla_CDS",
    "Aliame": "Aliame_CDS", "Codorniz": "Codorniz_CDS", "Gasacu": "Gasacu_CDS", "Masarm": "Masarm_CDS", "Perdiz": "Perdiz_CDS", "Tenrec": "Tenrec_CDS",
    "Alichi": "Alichi_CDS", "Colicorto": "Colicorto_CDS", "Gecko": "Gecko_CDS", "Melabo": "Melabo_CDS", "Perro": "Perro_CDS", "Toboba": "Toboba_CDS",
    "Ambrad": "Ambrad_CDS", "Colimbo": "Colimbo_CDS", "Geose": "Geose_CDS", "Microuni": "Microuni_CDS", "Plewal": "Plewal_CDS", "Torchina": "Torchina_CDS",
    "Avestruz": "Avestruz_CDS", "Conejo": "Conejo_CDS", "Grulla": "Grulla_CDS", "Mobhyp": "Mobhyp_CDS", "Podara": "Podara_CDS", "Tordia": "Tordia_CDS",
    "Boa": "Boa_CDS", "Cuco": "Cuco_CDS", "Hemicape": "Hemicape_CDS", "Monito": "Monito_CDS", "Porron": "Porron_CDS", "Torver": "Torver_CDS",
    "Bufalo": "Bufalo_CDS", "Dasnov": "Dasnov_CDS", "Hemoce": "Hemoce_CDS", "Mouse": "Mouse_CDS", "Pripec": "Pripec_CDS", "Tuatara": "Tuatara_CDS",
    "Bufo": "Bufo_CDS", "Demonio": "Demonio_CDS", "Humano": "Humano_CDS", "Mugcep": "Mugcep_CDS", "Proann": "Proann_CDS", "Tupaia": "Tupaia_CDS",
    "Bugar": "Bugar_CDS", "Dorexc": "Dorexc_CDS", "Hypsab": "Hypsab_CDS", "Murcielago": "Murcielago_CDS", "Punpun": "Punpun_CDS", "Vencejo": "Vencejo_CDS",
    "Caballo": "Caballo_CDS", "Elefante": "Elefante_CDS", "Inambu": "Inambu_CDS", "Musarana": "Musarana_CDS", "Ratem": "Ratem_CDS", "Xetro": "Xetro_CDS",
    "Calmil": "Calmil_CDS", "Elgaweb": "Elgaweb_CDS", "Kiwi": "Kiwi_CDS", "Nandu": "Nandu_CDS", "Rhina": "Rhina_CDS", "Zaraus": "Zaraus_CDS",
    "Carcar": "Carcar_CDS", "Emu": "Emu_CDS", "Latcha": "Latcha_CDS", "Ormela": "Ormela_CDS", "Rhityp": "Rhityp_CDS",
    "Chipla": "Chipla_CDS", "Equidna": "Equidna_CDS", "Lemur": "Lemur_CDS", "Ornitorrinco": "Ornitorrinco_CDS", "Scomax": "Scomax_CDS",
    "Chodid": "Chodid_CDS", "Funhet": "Funhet_CDS", "Leueri": "Leueri_CDS", "Patocollar": "PatoCollar_CDS", "Scycan": "Scycan_CDS",
    "Cisnene": "Cisnene_CDS", "Gallina": "Gallina_CDS", "Manati": "Manati_CDS", "Patozam": "Patozam_CDS", "Sorex": "Sorex_CDS"
}

# Extract ID from OGs files
def extraer_id(fasta_id):
    match = re.search(r'(XP|NP|YP|WP)_\d+\.\d+', fasta_id)  # Buscar XP_, NP_, etc.
    if match:
        return match.group(0)
    match = re.search(r'ENST\d+\.\d+', fasta_id) or re.search(r'ENSP\d+\.\d+', fasta_id)
    if match:
        return match.group(0)
    match = re.search(r'(\w+\.\d+)', fasta_id)  # Último recurso
    return match.group(0) if match else None

# Extract species names from OGs (they are in [])
def extraer_especie(fasta_id):
    matches = re.findall(r'\[(.*?)\]', fasta_id)
    return matches[-1] if matches else None

# Search the CDS sequence and the header
def buscar_cds(id_proteina, archivo_cds):
    if not os.path.exists(archivo_cds):
        print(f"⚠️ Archivo {archivo_cds} no encontrado.")
        return None, None
    
    with open(archivo_cds, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if id_proteina in record.description:
                return record.description, record.seq
    return None, None

# Print sequence and header if found
def procesar_fasta(archivo_proteinas):
    if not os.path.exists(archivo_proteinas):
        print(f"❌ Archivo {archivo_proteinas} no encontrado.")
        return
    
    with open(archivo_proteinas, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            id_proteina = extraer_id(record.description)
            especie = extraer_especie(record.description)
            
            if id_proteina and especie:
                especie_formateada = especie.replace(" ", "").capitalize()
                if especie_formateada in MAPA_ARCHIVOS:
                    archivo_cds = MAPA_ARCHIVOS[especie_formateada]
                    encabezado_cds, secuencia_cds = buscar_cds(id_proteina, archivo_cds)
                    
                    if encabezado_cds:
                        print(f">{encabezado_cds} [{especie}]")
                        print(secuencia_cds)
                    else:
                        print(f"❌ No se encontró {id_proteina} en {archivo_cds}")
                else:
                    print(f"⚠️ No hay archivo CDS asignado para la especie {especie} (formateado como {especie_formateada})")
            else:
                print(f"⚠️ No se pudo extraer ID o especie correctamente: {record.description}")

# Main arguments
def main():
    parser = argparse.ArgumentParser(description="Procesar un archivo FASTA de proteínas y buscar sus secuencias CDS.")
    parser.add_argument("archivo", help="Ruta al archivo FASTA de proteínas.")
    args = parser.parse_args()
    procesar_fasta(args.archivo)

if __name__ == "__main__":
    main()
