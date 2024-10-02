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

# Ribosome Binding Site (RBS) detection parameters
RBS_SEQUENCE = "AGGAGG"  # Shine-Dalgarno sequence
DEFAULT_RBS_UPSTREAM_RANGE = 20  # We will now check 4-20bp upstream of the start codon

def read_fna(filename):
    """Read the DNA sequence from a .fna file, removing header lines."""
    with open(filename, 'r') as file:
        lines = file.readlines()
    # Join lines together, excluding header (lines starting with '>')
    sequence = ''.join([line.strip() for line in lines if not line.startswith('>')])
    return sequence

def has_rbs(sequence, start_index, upstream_range=DEFAULT_RBS_UPSTREAM_RANGE):
    """Check for the presence of an RBS sequence upstream of a start codon."""
    upstream_region = sequence[max(0, start_index - 20):start_index - 4]
    return RBS_SEQUENCE in upstream_region

def translate_sequence(dna_sequence):
    """Translate a DNA sequence into an amino acid sequence using the codon table."""
    amino_acids = []
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i+3]
        amino_acids.append(codon_table.get(codon, ''))
    return ''.join(amino_acids)

def reverse_complement(sequence):
    """Get the reverse complement of a DNA sequence, including ambiguous bases."""
    complement = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K',
        'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', 'N': 'N'
    }
    return ''.join(complement.get(base, 'N') for base in reversed(sequence))

def find_genes_with_rbs(sequence, upstream_range=DEFAULT_RBS_UPSTREAM_RANGE, min_length=50):
    """Find ORFs in six reading frames with RBS filtering."""
    genes = []
    stop_codons = {'TAA', 'TAG', 'TGA'}
    
    # Check all three reading frames in the forward direction
    for frame in range(3):
        for i in range(frame, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if codon == 'ATG':  # Start codon
                if has_rbs(sequence, i, upstream_range=upstream_range):  # Check for RBS upstream
                    for j in range(i + 3, len(sequence) - 2, 3):
                        stop_codon = sequence[j:j+3]
                        if stop_codon in stop_codons:  # Check for stop codon
                            gene_seq = sequence[i:j+3]  # Get the gene sequence
                            if len(gene_seq) >= min_length * 3:  # Apply length filter (min_length codons)
                                amino_seq = translate_sequence(gene_seq)  # Translate DNA to amino acids
                                # Append the frame, start, stop, gene sequence, and amino sequence to the list
                                genes.append((frame + 1, i + 1, j + 3, gene_seq, amino_seq))
                            break  # Stop after finding the first valid gene in this frame
    
    # Also check all three reading frames in the reverse complement direction
    reverse_sequence = reverse_complement(sequence)
    for frame in range(3):
        for i in range(frame, len(reverse_sequence) - 2, 3):
            codon = reverse_sequence[i:i+3]
            if codon == 'ATG':  # Start codon
                if has_rbs(reverse_sequence, i, upstream_range=upstream_range):  # Check for RBS upstream
                    for j in range(i + 3, len(reverse_sequence) - 2, 3):
                        stop_codon = reverse_sequence[j:j+3]
                        if stop_codon in stop_codons:  # Check for stop codon
                            gene_seq = reverse_sequence[i:j+3]  # Get the gene sequence
                            if len(gene_seq) >= min_length * 3:  # Apply length filter (min_length codons)
                                amino_seq = translate_sequence(gene_seq)  # Translate DNA to amino acids
                                # Append the frame, start, stop, gene sequence, and amino sequence to the list
                                genes.append((-(frame + 1), i + 1, j + 3, gene_seq, amino_seq))
                            break  # Stop after finding the first valid gene in this frame
    
    return genes  # Return the list of identified genes

def main_task_6(filename, upstream_range=DEFAULT_RBS_UPSTREAM_RANGE):
    """Main function to read the sequence and find genes in six reading frames, filtering by RBS."""
    sequence = read_fna(filename)  # Read sequence from file
    genes = find_genes_with_rbs(sequence, upstream_range=upstream_range)  # Find genes with RBS
    
    # If genes were found, print their details
    if genes:
        print(f"Total number of ORFs found: {len(genes)}\n")
        for gene in genes:
            frame, start, stop, gene_seq, amino_seq = gene
            print(f"Frame {frame}: Start at {start}, Stop at {stop}")
            print(f"DNA Sequence: {gene_seq}")
            print(f"Amino Acid Sequence: {amino_seq}\n")
    else:
        print("No valid genes with RBS found.")  # Print message if no genes found

if __name__ == "__main__":
    fna_file = sys.argv[1]  # Take the filename from command-line arguments
    upstream_range = int(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_RBS_UPSTREAM_RANGE
    main_task_6(fna_file, upstream_range)  # Execute main function with the input file and optional upstream range
