import random
import argparse

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Append random integers to lines in a file.")
    parser.add_argument("--input_file", type=str, help="Path to the input samples file.")
    parser.add_argument("--output_file", type=str, help="Path to the output phenotype file.")
    parser.add_argument("--num_int", type=int, default=2, help="The number of categorical phenotypes to generate.")

    # Parse arguments
    args = parser.parse_args()

    file_in = args.input_file
    file_out = args.output_file
    num_int = args.num_int - 1

    with open(file_out, 'w') as g:
        with open(file_in, 'r') as f:
            for line in f:
                g.write(line.strip())
                g.write('\t')
                g.write(str(random.randint(0, num_int)))
                g.write('\n')

if __name__ == "__main__":
    main()
