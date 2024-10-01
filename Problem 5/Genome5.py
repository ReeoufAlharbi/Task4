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
    'TAC':'Y', 'TAT':'Y',                        # Tyrosine
    'TGC':'C', 'TGT':'C',                        # Cysteine
    'TGG':'W',                                   # Tryptophan
    'TAA':'', 'TAG':'', 'TGA':''                 # Stop codons (empty strings as placeholders)
}

def read_fna(filename):
    """Read the DNA sequence from a .fna file, removing header lines."""
    with open(filename, 'r') as file:
        lines = file.readlines()
    # Join lines together, excluding header (lines starting with '>')
    sequence = ''.join([line.strip() for line in lines if not line.startswith('>')])
    return sequence

def reverse_complement(sequence):
    """Generate the reverse complement of a DNA sequence, handling unknown nucleotides."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}  # Handle 'N' by mapping it to 'N'
    return ''.join(complement.get(base, 'N') for base in reversed(sequence))  # Use 'N' as a fallback for unknown bases

def translate_sequence(dna_seq):
    """Translate a DNA sequence into an amino acid sequence."""
    amino_acids = []
    for i in range(0, len(dna_seq) - 2, 3):  # Process in triplets (codons)
        codon = dna_seq[i:i+3]
        amino_acids.append(codon_table.get(codon, '?'))  # Translate codon, '?' for unknown codons
    return ''.join(amino_acids)

def find_genes_six_frames(sequence, min_length=100):
    """Find genes in all six reading frames with a filter for minimum length."""
    genes = []
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]

    # Search in all three forward frames
    for frame in range(3):
        for i in range(frame, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if codon == start_codon:
                for j in range(i + 3, len(sequence) - 2, 3):
                    stop_codon = sequence[j:j+3]
                    if stop_codon in stop_codons:
                        gene_seq = sequence[i:j+3]  # Extract gene sequence
                        if len(gene_seq) >= min_length * 3:  # Apply length filter (min_length codons)
                            amino_seq = translate_sequence(gene_seq)  # Translate gene sequence
                            genes.append(('Forward', frame + 1, i + 1, j + 3, gene_seq, amino_seq))
                        break  # Stop after finding the first valid gene in this frame

    reverse_seq = reverse_complement(sequence)  # Get reverse complement of the sequence

    # Search in all three reverse frames
    for frame in range(3):
        for i in range(frame, len(reverse_seq) - 2, 3):
            codon = reverse_seq[i:i+3]
            if codon == start_codon:
                for j in range(i + 3, len(reverse_seq) - 2, 3):
                    stop_codon = reverse_seq[j:j+3]
                    if stop_codon in stop_codons:
                        start_pos = len(sequence) - (i + 1)  # Convert reverse position to original sequence position
                        stop_pos = len(sequence) - (j + 3)
                        gene_seq = reverse_seq[i:j+3]  # Extract gene sequence
                        if len(gene_seq) >= min_length * 3:  # Apply length filter (min_length codons)
                            amino_seq = translate_sequence(gene_seq)  # Translate gene sequence
                            genes.append(('Reverse', frame + 1, stop_pos + 1, start_pos + 1, gene_seq, amino_seq))
                        break  # Stop after finding the first valid gene in this frame

    return genes

def main_task_5(filename):
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
        print("No valid start-stop codon region found.")

if __name__ == "__main__":
    fna_file = sys.argv[1]  # Take the filename from command-line arguments
    main_task_5(fna_file)  # Execute main function with the input file
