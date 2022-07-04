import argparse

class ParsedInputArguments():
    def __init__(self) -> None:
        parser = argparse.ArgumentParser(
            description='Targeton Designer CLI')

        parser.add_argument('command',
                            help='Command to run in Designer CLI',
                            type=str)

        parser = add_input_args(parser)

        self.arguments = parser.parse_args()
        self.command = self.arguments.command

    def get_command(self) -> str:
        return self.command

    def get_args(self) -> str:
        return self.arguments

def positive_int(arg: str) -> int:
    if int(arg) <= 0:
        raise argparse.ArgumentTypeError('Parameter must be above 0')
    return int(arg)

def len_positive_int(arg: str) -> int:
    if 10000 < int(arg) or int(arg) <= 0:
        raise argparse.ArgumentTypeError('Parameter must be above 0 and below 10000')
    return int(arg)

def add_input_args(parser):
    parser.add_argument('--bed',
                        help='BED file from the slicer tool containing coords, strand and IDs')
    parser.add_argument('--fasta',
                        help='FASTA file from the slicer tool containing seqs and IDs')
    parser.add_argument('-f5', '--flank_5',
                        help='how far to extend region at 5\' end (default 50nt)',
                        type=int, default=50)
    parser.add_argument('-f3', '--flank_3',
                        help='how far to extend region at 3\' end (default 50nt)',
                        type=int, default=50)
    parser.add_argument('-l', '--length',
                        help='length of each slice (default 210nt)',
                        type=len_positive_int, default=210)
    parser.add_argument('-o', '--offset',
                        help='offset between each slice (default 5nt)',
                        type=positive_int, default=5)
    parser.add_argument('-d', '--dir',
                        help='output directory name to be timestamped'
                             ' (default \'td_output\')',
                        type=str, default='td_output')
    
    return parser

    