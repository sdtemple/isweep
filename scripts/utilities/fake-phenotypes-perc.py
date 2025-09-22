import random
import argparse

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Append random integers to lines in a file.")
    parser.add_argument("--input_file", type=str, help="Path to the input samples file.")
    parser.add_argument("--output_file", type=str, help="Path to the output phenotype file.")
    parser.add_argument("--percent_base", type=float, default=0.5, help="Percentage of base category.")

    # Parse arguments
    args = parser.parse_args()

    file_in = args.input_file
    file_out = args.output_file
    percent_base = args.percent_base

    with open(file_out, 'w') as g:
        with open(file_in, 'r') as f:
            for line in f:
                g.write(line.strip())
                g.write('\t')
                randvalue = random.random()
                if randvalue <= percent_base:
                    g.write('1')
                else:
                    g.write('0')
                g.write('\n')

if __name__ == "__main__":
    main()
