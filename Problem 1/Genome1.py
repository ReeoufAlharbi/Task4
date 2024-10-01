import sys  # Importing the sys module to handle command-line arguments

# Function to read a fna file
def read_fna(filename):
    # Open the given file in read mode
    with open(filename, 'r') as file:
        # Read all lines from the file into a list
        lines = file.readlines()
    # Join the sequence lines while ignoring any lines that start with '>' (fna headers)
    sequence = ''.join([line.strip() for line in lines if not line.startswith('>')])
    return sequence  # Return the joined sequence without headers or line breaks

# Function to find genes in three reading frames
def find_genes(sequence):
    start_codon = 'ATG'  # Define the start codon for gene translation (ATG)
    stop_codons = ['TAA', 'TAG', 'TGA']  # Define the three possible stop codons
    genes = []  # Initialize an empty list to store the identified genes
    
    # Loop through the three reading frames (starting at frame 0, 1, and 2)
    for frame in range(3):
        # Loop through the sequence in steps of 3 nucleotides (codons), starting at the current frame
        for i in range(frame, len(sequence) - 2, 3):
            codon = sequence[i:i+3]  # Extract a codon (3 nucleotides)
            if codon == start_codon:  # Check if the codon is the start codon (ATG)
                # Once a start codon is found, look for the next stop codon in the same reading frame
                for j in range(i + 3, len(sequence) - 2, 3):
                    stop_codon = sequence[j:j+3]  # Extract another codon (stop codon check)
                    if stop_codon in stop_codons:  # If a stop codon is found
                        # Store the frame (1, 2, or 3), the start and stop positions, and the gene sequence
                        genes.append((frame + 1, i + 1, j + 3, sequence[i:j+3]))
                        break  # Break the inner loop after finding a valid gene
    return genes  # Return the list of identified genes

# Main function for Task 1
def main_task_1(filename):
    sequence = read_fna(filename)  # Read the sequence from the fna file
    genes = find_genes(sequence)  # Find genes in the sequence
    if genes:  # If any genes are found
        # Loop through the identified genes and print the reading frame, start/stop positions, and gene sequence
        for gene in genes:
            print(f"Frame {gene[0]}: Start at {gene[1]}, Stop at {gene[2]}: {gene[3]}")
    else:
        # If no genes are found, print an informative message
        print("No valid start-stop codon region found.")

# Usage example: python script.py input.fna
if __name__ == "__main__":
    fna_file = sys.argv[1]  # Take the fna file as a command-line argument
    main_task_1(fna_file)  # Run the main function using the input fna file
