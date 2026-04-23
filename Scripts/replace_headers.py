# Original Fasta
fasta_file = "../Path_to_proteomes/Mouse.faa"

# File with new IDs
nuevas_ids_file = "Mouse_ids_4"

# Output
output_file = "Mouse.fa"

# Read new IDs from file
nuevas_ids = []
with open(nuevas_ids_file, "r") as ids_file:
    nuevas_ids = [line.strip() for line in ids_file if line.strip()]

# Process Fasta file and the new IDs file
with open(fasta_file, "r") as fasta, open(output_file, "w") as salida:
    nueva_id_actual = 0
    for line in fasta:
        if line.startswith(">"):
            # Replace the ID in the actual line with the new ID
            nueva_id = nuevas_ids[nueva_id_actual]
            salida.write(f">{nueva_id}\n")
            nueva_id_actual += 1
        else:
            # Keep the sequence
            salida.write(line)

# Print that the process is finished
print("Reemplazo de IDs completado. El archivo modificado se ha guardado en", output_file)
