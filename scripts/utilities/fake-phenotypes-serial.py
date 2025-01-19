import random
import argparse

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Append serially integers to lines in a file.")
    parser.add_argument("--input_file", type=str, help="Path to the input samples file.")
    parser.add_argument("--output_file", type=str, help="Path to the output phenotype file.")
    parser.add_argument("--num_first_serial", type=int, default=2, help="First this many lines get 1. Rest get 0.")

    # Parse arguments
    args = parser.parse_args()

    file_in = args.input_file
    file_out = args.output_file
    num_int = args.num_first_serial

    ctr = 0
    with open(file_out, 'w') as g:
        with open(file_in, 'r') as f:
            for line in f:
                g.write(line.strip())
                g.write('\t')
                if ctr <= num_int:
                    g.write('1\t')
                else:
                    g.write('0\t')
                g.write('\n')

if __name__ == "__main__":
    main()
