import sys  # Importing sys to handle command-line arguments

# Codon to amino acid translation table
# This dictionary maps each codon (3-letter DNA sequence) to its corresponding amino acid.
codon_table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',  # Isoleucine and Methionine
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',  # Threonine
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',  # Asparagine, Lysine
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',  # Serine, Arginine
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',  # Leucine
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',  # Proline
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',  # Histidine, Glutamine
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',  # Arginine
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',  # Valine
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',  # Alanine
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',  # Aspartic acid, Glutamic acid
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',  # Glycine
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',  # Serine
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',  # Phenylalanine, Leucine
    'TAC':'Y', 'TAT':'Y',                        # Tyrosine
    'TGC':'C', 'TGT':'C',                        # Cysteine
    'TGG':'W',                                   # Tryptophan
    'TAA':'', 'TAG':'', 'TGA':''                 # Stop codons (empty strings as placeholders)
}

def read_fna(filename):
    """Reads a .fna file and extracts the DNA sequence."""
    with open(filename, 'r') as file:  # Open the file for reading
        lines = file.readlines()       # Read all lines of the file
    # Join all lines that don't start with '>' (which represents metadata) to form the sequence
    sequence = ''.join([line.strip() for line in lines if not line.startswith('>')])
    return sequence  # Return the DNA sequence

def reverse_complement(sequence):
    """Generates the reverse complement of a given DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}  # Base pairing rules
    # Reverse the sequence and replace each base with its complement
    return ''.join(complement[base] for base in reversed(sequence))

def translate_sequence(dna_seq):
    """Translates a DNA sequence into an amino acid sequence."""
    amino_acids = []  # List to store amino acids
    # Loop through the DNA sequence in chunks of 3 (codons)
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]  # Extract a codon
        amino_acids.append(codon_table.get(codon, '?'))  # Translate codon to amino acid or '?' if unknown
    return ''.join(amino_acids)  # Return the translated amino acid sequence

def find_genes_six_frames(sequence):
    """Finds open reading frames (ORFs) in all 6 frames and translates them into amino acid sequences."""
    start_codon = 'ATG'  # Start codon for translation initiation
    stop_codons = ['TAA', 'TAG', 'TGA']  # Stop codons for translation termination
    genes = set()  # Use a set to store amino acid sequences (to remove duplicates)
    
    # Forward direction (3 reading frames)
    for frame in range(3):  # Shift by 0, 1, or 2 bases
        for i in range(frame, len(sequence) - 2, 3):  # Read codons in the current frame
            codon = sequence[i:i+3]
            if codon == start_codon:  # Found a start codon
                # Search for the next stop codon
                for j in range(i + 3, len(sequence) - 2, 3):
                    stop_codon = sequence[j:j+3]
                    if stop_codon in stop_codons:  # Found a stop codon
                        gene_seq = sequence[i:j+3]  # Extract the gene sequence
                        amino_seq = translate_sequence(gene_seq)  # Translate the gene to amino acids
                        genes.add(amino_seq)  # Add the amino acid sequence to the set
                        break  # Stop searching once a stop codon is found

    # Reverse direction (3 reading frames on the reverse complement)
    reverse_seq = reverse_complement(sequence)  # Generate the reverse complement of the DNA sequence
    for frame in range(3):  # Shift by 0, 1, or 2 bases
        for i in range(frame, len(reverse_seq) - 2, 3):  # Read codons in the current frame
            codon = reverse_seq[i:i+3]
            if codon == start_codon:  # Found a start codon
                # Search for the next stop codon
                for j in range(i + 3, len(reverse_seq) - 2, 3):
                    stop_codon = reverse_seq[j:j+3]
                    if stop_codon in stop_codons:  # Found a stop codon
                        gene_seq = reverse_seq[i:j+3]  # Extract the gene sequence
                        amino_seq = translate_sequence(gene_seq)  # Translate the gene to amino acids
                        genes.add(amino_seq)  # Add the amino acid sequence to the set
                        break  # Stop searching once a stop codon is found
                        
    return genes  # Return the set of unique amino acid sequences

def main_task_3(filename):
    """Reads a DNA sequence, finds all ORFs, and prints unique amino acid sequences."""
    sequence = read_fna(filename)  # Read the DNA sequence from the .fna file
    amino_acid_sequences = find_genes_six_frames(sequence)  # Find ORFs in all 6 frames
    
    if amino_acid_sequences:
        # Print the unique amino acid sequences (ORFs found)
        for amino_seq in amino_acid_sequences:
            print(f"{amino_seq}")
    else:
        print("No valid start-stop codon region found.")  # Message if no ORFs are found

# Main program entry point
if __name__ == "__main__":
    fna_file = sys.argv[1]  # Take the filename as the first command-line argument
    main_task_3(fna_file)  # Run the task with the provided file
