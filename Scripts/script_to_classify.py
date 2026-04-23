#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import sys
import os
from collections import defaultdict

# Dictionary for classification species

clasificacion = {
    "Humano": {"clase": "Mamíferos", "grupo": "Euarchontoglires"},
    "Mouse": {"clase": "Mamíferos", "grupo": "Euarchontoglires"},
    "Lemur": {"clase": "Mamíferos", "grupo": "Euarchontoglires"},
    "Tupaia": {"clase": "Mamíferos", "grupo": "Euarchontoglires"},
    "Conejo": {"clase": "Mamíferos", "grupo": "Euarchontoglires"},
    "Perro": {"clase": "Mamíferos", "grupo": "Laurasiatheria"},
    "Murcielago": {"clase": "Mamíferos", "grupo": "Laurasiatheria"},
    "Sorex": {"clase": "Mamíferos", "grupo": "Laurasiatheria"},
    "Caballo": {"clase": "Mamíferos", "grupo": "Laurasiatheria"},
    "Bufalo": {"clase": "Mamíferos", "grupo": "Laurasiatheria"},
    "Elefante": {"clase": "Mamíferos", "grupo": "Atlantogenata"},
    "Manati": {"clase": "Mamíferos", "grupo": "Atlantogenata"},
    "Tenrec": {"clase": "Mamíferos", "grupo": "Atlantogenata"},
    "Dasnov": {"clase": "Mamíferos", "grupo": "Atlantogenata"},
    "Musarana": {"clase": "Mamíferos", "grupo": "Atlantogenata"},
    "Chodid": {"clase": "Mamíferos", "grupo": "Atlantogenata"},
    "Colicorto": {"clase": "Mamíferos", "grupo": "Marsupialia"},
    "Demonio": {"clase": "Mamíferos", "grupo": "Marsupialia"},
    "Monito": {"clase": "Mamíferos", "grupo": "Marsupialia"},
    "Zaraus": {"clase": "Mamíferos", "grupo": "Marsupialia"},
    "Mardo": {"clase": "Mamíferos", "grupo": "Marsupialia"},
    "Equidna": {"clase": "Mamíferos", "grupo": "Monotrema"},
    "Ornitorrinco": {"clase": "Mamíferos", "grupo": "Monotrema"},
    "Gallo": {"clase": "Aves", "grupo": "Galliformes"},
    "Codorniz": {"clase": "Aves", "grupo": "Galliformes"},
    "Pavo": {"clase": "Aves", "grupo": "Galliformes"},
    "Gallina": {"clase": "Aves", "grupo": "Galliformes"},
    "Perdiz": {"clase": "Aves", "grupo": "Galliformes"},
    "PatoCollar": {"clase": "Aves", "grupo": "Anseriformes"},
    "Cisnevul": {"clase": "Aves", "grupo": "Anseriformes"},
    "Cisnene": {"clase": "Aves", "grupo": "Anseriformes"},
    "Porron": {"clase": "Aves", "grupo": "Anseriformes"},
    "Patozam": {"clase": "Aves", "grupo": "Anseriformes"},
    "Inambu": {"clase": "Aves", "grupo": "Paleognathos"},
    "Kiwi": {"clase": "Aves", "grupo": "Paleognathos"},
    "Emu": {"clase": "Aves", "grupo": "Paleognathos"},
    "Nandu": {"clase": "Aves", "grupo": "Paleognathos"},
    "Avestruz": {"clase": "Aves", "grupo": "Paleognathos"},
    "Vencejo": {"clase": "Aves", "grupo": "Neoaves"},
    "Cuco": {"clase": "Aves", "grupo": "Neoaves"},
    "Grulla": {"clase": "Aves", "grupo": "Neoaves"},
    "Colimbo": {"clase": "Aves", "grupo": "Neoaves"},
    "Aguila": {"clase": "Aves", "grupo": "Neoaves"},
    "Tordia": {"clase": "Reptiles", "grupo": "Testudines"},
    "Toboba": {"clase": "Reptiles", "grupo": "Testudines"},
    "Torchina": {"clase": "Reptiles", "grupo": "Testudines"},
    "Torver": {"clase": "Reptiles", "grupo": "Testudines"},
    "Aliame": {"clase": "Reptiles", "grupo": "Crocodilia"},
    "Alichi": {"clase": "Reptiles", "grupo": "Crocodilia"},
    "Hemicape": {"clase": "Reptiles", "grupo": "Squamata"},
    "Elgaweb": {"clase": "Reptiles", "grupo": "Squamata"},
    "Gecko": {"clase": "Reptiles", "grupo": "Squamata"},
    "Podara": {"clase": "Reptiles", "grupo": "Squamata"},
    "Boa": {"clase": "Reptiles", "grupo": "Squamata"},
    "Tuatara": {"clase": "Reptiles", "grupo": "Rhynchocephalia"},
    "Xetro": {"clase": "Anfibios", "grupo": "Anura"},
    "Ratem": {"clase": "Anfibios", "grupo": "Anura"},
    "Bugar": {"clase": "Anfibios", "grupo": "Anura"},
    "Bufo": {"clase": "Anfibios", "grupo": "Anura"},
    "Plewal": {"clase": "Anfibios", "grupo": "Caudata"},
    "Geose": {"clase": "Anfibios", "grupo": "Gymnophiona"},
    "Microuni": {"clase": "Anfibios", "grupo": "Gymnophiona"},
    "Rhina": {"clase": "Anfibios", "grupo": "Gymnophiona"},
    "Melabo": {"clase": "Oseos", "grupo": "Actinopterygii"},
    "Ormela": {"clase": "Oseos", "grupo": "Actinopterygii"},
    "Funhet": {"clase": "Oseos", "grupo": "Actinopterygii"},
    "Gasacu": {"clase": "Oseos", "grupo": "Actinopterygii"},
    "Mugcep": {"clase": "Oseos", "grupo": "Actinopterygii"},
    "Punpun": {"clase": "Oseos", "grupo": "Actinopterygii"},
    "Scomax": {"clase": "Oseos", "grupo": "Actinopterygii"},
    "Masarm": {"clase": "Oseos", "grupo": "Actinopterygii"},
    "Dorexc": {"clase": "Oseos", "grupo": "Actinopterygii"},
    "Takfla": {"clase": "Oseos", "grupo": "Actinopterygii"},
    "Latcha": {"clase": "Oseos", "grupo": "Actinistia"},
    "Proann": {"clase": "Oseos", "grupo": "Dipnomorpha"},
    "Leueri": {"clase": "Cartilaginosos", "grupo": "Batoidea"},
    "Pripec": {"clase": "Cartilaginosos", "grupo": "Batoidea"},
    "Mobhyp": {"clase": "Cartilaginosos", "grupo": "Batoidea"},
    "Hypsab": {"clase": "Cartilaginosos", "grupo": "Batoidea"},
    "Ambrad": {"clase": "Cartilaginosos", "grupo": "Batoidea"},
    "Carcar": {"clase": "Cartilaginosos", "grupo": "Selachii"},
    "Scycan": {"clase": "Cartilaginosos", "grupo": "Selachii"},
    "Chipla": {"clase": "Cartilaginosos", "grupo": "Selachii"},
    "Rhityp": {"clase": "Cartilaginosos", "grupo": "Selachii"},
    "Hemoce": {"clase": "Cartilaginosos", "grupo": "Selachii"},
    "Calmil": {"clase": "Cartilaginosos", "grupo": "Chimaeriformes"},
    
}


# Orden global predefinido
orden_global = []
for clase in sorted(set(d["clase"] for d in clasificacion.values())):
    grupos_clase = sorted(set(g["grupo"] for g in clasificacion.values() if g["clase"] == clase))
    orden_global.extend(grupos_clase)
orden_global = list(dict.fromkeys(orden_global))  # Deduplicar

def procesar_archivo(archivo_entrada):
    """
    Process input file and count species per group

    Args:
        Input file (str): Path to input file.

    Returns:
        dict: Count.
    """
    conteo = defaultdict(int)
    patron = re.compile(r'\[([^\[\]]+)\]\s*$')

    try:
        with open(archivo_entrada, 'r', encoding='utf-8') as f:
            for linea in f:
                linea = linea.strip()
                if linea.startswith(">"):  # IDs
                    match = patron.search(linea)
                    if match:
                        especie = match.group(1).strip()
                        if especie in clasificacion:
                            grupo = clasificacion[especie]["grupo"]
                            conteo[grupo] += 1
        return conteo
    except FileNotFoundError:
        print(f"El archivo '{archivo_entrada}' no fue encontrado.")
        return {}
    except Exception as e:
        print(f"Error al procesar '{archivo_entrada}': {e}")
        return {}

def imprimir_resultados_con_formato(conteo, imprimir_encabezados=True):
    """
    Print results with the specified format.

    Args:
        conteo (dict): Dictionary with the count of species for group.
        imprimir_encabezados (bool): if you want to print the header.
    """
    if imprimir_encabezados:
        print(";".join(orden_global))  # Print headers

    # Generate the counts in order
    conteos_ordenados = [str(conteo.get(grupo, 0)) for grupo in orden_global]
    print(";".join(conteos_ordenados))

def main():
    """
    Main function.
    """
    # Read the input file with the specified OMAGroups to use
    lista_archivos = "OG_List"
    if not os.path.exists(lista_archivos):
        print(f"No se encontró el archivo '{lista_archivos}'.")
        return

    with open(lista_archivos, 'r', encoding='utf-8') as f:
        archivos_a_procesar = [linea.strip() for linea in f if linea.strip()]

    # Process each file in our list
    for idx, archivo in enumerate(archivos_a_procesar):
        if not os.path.exists(archivo):
            print(f"El archivo '{archivo}' no fue encontrado, se omitirá.")
            continue

        conteo = procesar_archivo(archivo)

        # Print the header
        imprimir_encabezados = idx == 0
        imprimir_resultados_con_formato(conteo, imprimir_encabezados=imprimir_encabezados)

if __name__ == "__main__":
    main()
