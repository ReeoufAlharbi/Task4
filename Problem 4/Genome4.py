import sys  # Import sys to handle command-line arguments

# Codon to amino acid translation table
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
    'TAC':'Y', 'TAT':'Y',  # Tyrosine
    'TGC':'C', 'TGT':'C',  # Cysteine
    'TGG':'W',  # Tryptophan
    'TAA':'', 'TAG':'', 'TGA':''  # Stop codons (empty strings as placeholders)
}

def read_fna(filename):
    """Read the DNA sequence from a .fna file, removing header lines."""
    with open(filename, 'r') as file:  # Open the file for reading
        lines = file.readlines()  # Read all lines in the file
    # Join lines together, excluding header (lines starting with '>')
    sequence = ''.join([line.strip() for line in lines if not line.startswith('>')])
    return sequence  # Return the cleaned-up DNA sequence

def reverse_complement(sequence):
    """Generate the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}  # Define nucleotide complements
    # Create reverse complement by reversing and substituting bases
    return ''.join([complement[base] for base in reversed(sequence)])

def translate_sequence(sequence):
    """Translate a DNA sequence into an amino acid sequence using the codon table."""
    amino_acids = []
    # Iterate over the sequence in steps of 3 (codon length)
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]  # Extract codon
        amino_acids.append(codon_table.get(codon, ''))  # Translate codon into amino acid, default to '' for stop codons
    return ''.join(amino_acids)  # Return the amino acid sequence

def find_genes_six_frames(sequence):
    """Find genes in all six reading frames (three forward, three reverse)."""
    start_codon = 'ATG'  # Start codon is ATG
    stop_codons = {'TAA', 'TAG', 'TGA'}  # Define stop codons
    genes = []  # List to store gene information

    # Search in all three forward frames
    for frame in range(3):
        for i in range(frame, len(sequence) - 2, 3):  # Iterate over sequence in steps of 3
            codon = sequence[i:i+3]  # Extract codon
            if codon == start_codon:  # Check if it's a start codon
                for j in range(i + 3, len(sequence) - 2, 3):  # Search for stop codon after start
                    stop_codon = sequence[j:j+3]  # Extract stop codon
                    if stop_codon in stop_codons:  # If stop codon is found
                        gene_seq = sequence[i:j+3]  # Extract the gene sequence
                        amino_seq = translate_sequence(gene_seq)  # Translate gene to amino acid sequence
                        # Store gene information: strand, frame, start, stop, gene sequence, amino acid sequence
                        genes.append(('Forward', frame + 1, i + 1, j + 3, gene_seq, amino_seq))
                        break  # Stop after finding the first valid gene in this frame

    reverse_seq = reverse_complement(sequence)  # Get reverse complement of the sequence

    # Search in all three reverse frames
    for frame in range(3):
        for i in range(frame, len(reverse_seq) - 2, 3):  # Iterate over reverse sequence
            codon = reverse_seq[i:i+3]  # Extract codon
            if codon == start_codon:  # Check for start codon in reverse
                for j in range(i + 3, len(reverse_seq) - 2, 3):  # Search for stop codon
                    stop_codon = reverse_seq[j:j+3]  # Extract stop codon
                    if stop_codon in stop_codons:  # If stop codon is found
                        start_pos = len(sequence) - (i + 1)  # Convert reverse position to original sequence position
                        stop_pos = len(sequence) - (j + 3)
                        gene_seq = reverse_seq[i:j+3]  # Extract gene sequence
                        amino_seq = translate_sequence(gene_seq)  # Translate gene sequence
                        # Store gene information for reverse strand
                        genes.append(('Reverse', frame + 1, stop_pos + 1, start_pos + 1, gene_seq, amino_seq))
                        break  # Stop after finding the first valid gene in this frame

    return genes  # Return list of all found genes

def main_task_4(filename):
    """Main function to read the sequence and find genes in six reading frames."""
    sequence = read_fna(filename)  # Read sequence from file
    genes = find_genes_six_frames(sequence)  # Find genes in all six frames
    if genes:
        for gene in genes:
            strand, frame, start, stop, gene_seq, amino_seq = gene
            # Output the results with gene details
            print(f"{strand} Strand Frame {frame}: Start at {start}, Stop at {stop}")
            print(f"DNA Sequence: {gene_seq}")
            print(f"Amino Acid Sequence: {amino_seq}\n")
    else:
        print("No valid start-stop codon region found.")  # Message if no genes are found

if __name__ == "__main__":
    fna_file = sys.argv[1]  # Take the filename from command-line arguments
    main_task_4(fna_file)  # Execute main function with the input file
