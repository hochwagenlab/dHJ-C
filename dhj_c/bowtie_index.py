import optparse
import subprocess

def build_bowtie_index(fa_input, output_prefix):
    """
    Build Bowtie1 index at output_prefix using fa_input
    fa_input can be single fasta file or a comma delimited list of files
    """
    subprocess.check_call([
        'bowtie-build',
        fa_input,
        output_prefix,
    ])

if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('-f', dest='fa_input')
    parser.add_option('-o', dest='output_prefix')
    options, args = parser.parse_args()
    build_bowtie_index(
        fa_input=options.fa_input, 
        output_prefix=options.output_prefix,
    )

