import sys

# Function to read a .fna file (FASTA format) and return the sequence without header lines
def read_fna(filename):
    with open(filename, 'r') as file:  # Open the file for reading
        lines = file.readlines()  # Read all lines from the file
    # Join all lines, excluding those that start with '>', and remove whitespace characters
    sequence = ''.join([line.strip() for line in lines if not line.startswith('>')])
    return sequence  # Return the cleaned DNA sequence

# Function to compute the reverse complement of a DNA sequence
def reverse_complement(sequence):
    # Dictionary to map each nucleotide to its complement
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    # Reverse the sequence and replace each base with its complement
    return ''.join(complement[base] for base in reversed(sequence))

# Function to find genes in all six reading frames (3 forward, 3 reverse)
def find_genes_six_frames(sequence):
    start_codon = 'ATG'  # The standard start codon for protein synthesis
    stop_codons = ['TAA', 'TAG', 'TGA']  # Standard stop codons in DNA sequences
    genes = []  # Initialize a list to store found genes
    
    # Loop over the three possible forward reading frames
    for frame in range(3):
        # Loop through the sequence in steps of 3 (codon length)
        for i in range(frame, len(sequence) - 2, 3):
            codon = sequence[i:i+3]  # Extract a codon (3 nucleotides)
            if codon == start_codon:  # If start codon is found
                # Search for the next stop codon
                for j in range(i + 3, len(sequence) - 2, 3):
                    stop_codon = sequence[j:j+3]
                    if stop_codon in stop_codons:  # If a stop codon is found
                        # Record the gene details: strand, frame, start, stop, and gene sequence
                        genes.append(('Forward', frame + 1, i + 1, j + 3, sequence[i:j+3]))
                        break  # Stop looking after finding the first valid stop codon

    # Compute the reverse complement of the sequence for reverse frames
    reverse_seq = reverse_complement(sequence)

    # Loop over the three possible reverse reading frames
    for frame in range(3):
        for i in range(frame, len(reverse_seq) - 2, 3):
            codon = reverse_seq[i:i+3]
            if codon == start_codon:
                # Search for the next stop codon
                for j in range(i + 3, len(reverse_seq) - 2, 3):
                    stop_codon = reverse_seq[j:j+3]
                    if stop_codon in stop_codons:
                        # Calculate the position in the original sequence
                        start_pos = len(sequence) - (i + 1)
                        stop_pos = len(sequence) - (j + 3)
                        # Record the gene details for the reverse strand
                        genes.append(('Reverse', frame + 1, stop_pos + 1, start_pos + 1, reverse_seq[i:j+3]))
                        break  # Stop after finding the first valid stop codon
    
    return genes  # Return the list of found genes

# Main function to handle input and output
def main_task_2(filename):
    sequence = read_fna(filename)  # Read the DNA sequence from the .fna file
    genes = find_genes_six_frames(sequence)  # Find all genes in the sequence
    if genes:  # If any genes were found
        for gene in genes:
            # Unpack the gene details: strand, frame, start, stop, and gene sequence
            strand, frame, start, stop, gene_seq = gene
            # Print the details of each found gene
            print(f"{strand} Strand Frame {frame}: Start at {start}, Stop at {stop}: {gene_seq}")
    else:
        # If no genes were found, print a message indicating that
        print("No valid start-stop codon region found.")

# Entry point of the script
if __name__ == "__main__":
    fna_file = sys.argv[1]  # Get the input filename from command line arguments
    main_task_2(fna_file)  # Call the main function with the filename
