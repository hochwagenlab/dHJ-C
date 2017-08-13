import optparse
import subprocess
import tempfile

from Bio import SeqIO


def align(fastq_in, ):
    """
    Align fastq file using Bowtie
    Remove unmapped reads using Samtools
    """
    subprocess.check_call([
        'bowtie',
        '-q',
        '-m', '1',
        '-v', '0',
        '-p', '8',
        '-S',  # SAM output
        '-3', '129', '-5', '2',  # Read Trimming
        '--un', unaligned_fastq,
        '--max', multi_fastq,
        bowtie_index,
        '-a', input_fastq,
        sam_output
    ])
    subprocess.check_call([
        'samtools', 'view',
        '-F', '4',
        '-h',
        '-S',
        '-o', 'Bowtie/${TAG1}_SK1-${i}_R1-Align1.sam',
        'Bowtie/${TAG1}_SK1-${i}_R1.sam'
    ])


def read_merge(lib_name, fastq_num, read):
    """Merge alignment lists from SAM files, covering combinations of
    different genomic templates for R1 and R2."""
    # Read in S288c alignments
    with open(lib_name + '_S288c-' + fastq_num + '_' + read + \
              '-data5.txt') as file:
        contents = file.readlines()
    # make S288c chromosome notation compatible
    chr = {'chrI':1, 'chrII':2, 'chrIII':3, 'chrIV':4, 'chrV':5,
                  'chrVI':6, 'chrVII':7, 'chrVIII':8, 'chrIX':9, 'chrX':10,
                  'chrXI':11, 'chrXII':12, 'chrXIII':13, 'chrXIV':14,
                  'chrXV':15, 'chrXVI':16}
    # split lines and make dictionary of unique reads
    unique_reads = {}
    for line in contents:
        tabbed = line.replace('_', '\t')
        aln = tabbed.strip().split('\t')
        unique_reads[aln[0]] = [aln[1], str(chr[aln[2]]), aln[3], aln[5]]
    print str(len(unique_reads.keys())) + ' ' + read + ' reads align to S288c'
    # repeat for SK1
    # Read in SK1 alignments
    with open(lib_name + '_SK1-' + fastq_num + '_' + read + \
              '-data5.txt') as file:
        contents = file.readlines()  
    # make SK1 chromosome notation compatible
    chr = {'chr01':1, 'chr02':2, 'chr03':3, 'chr04':4, 'chr05':5,
                   'chr06':6, 'chr07':7, 'chr08':8, 'chr09':9, 'chr10':10,
                   'chr11':11, 'chr12':12, 'chr13':13, 'chr14':14,
                   'chr15':15, 'chr16':16}
    # split lines and add to dictionary of unique reads
    for line in contents:
        tabbed = line.replace('_', '\t')
        aln = tabbed.strip().split('\t')
        if unique_reads.get(aln[0]):
            continue
        else:
            unique_reads[aln[0]] = [aln[1], str(chr[aln[2]]), aln[3], aln[5]]
    print str(len(unique_reads.keys())) + ' ' + read + \
               ' reads align to S288c or SK1'
    return unique_reads

def read_merge_slow(lib_name, fastq_num, read):
    """Merge alignment lists from SAM files, covering combinations of
    different genomic templates for R1 and R2."""
    # Read in S288c alignments
    # make S288c chromosome notation compatible
    chr = {'chrI':1, 'chrII':2, 'chrIII':3, 'chrIV':4, 'chrV':5,
                  'chrVI':6, 'chrVII':7, 'chrVIII':8, 'chrIX':9, 'chrX':10,
                  'chrXI':11, 'chrXII':12, 'chrXIII':13, 'chrXIV':14,
                  'chrXV':15, 'chrXVI':16}
    # Read in S288c alignments, line by line
    f = open(lib_name + '_S288c-' + fastq_num + '_' + read + '-data5.txt',
             'r')
    unique_reads = {}
    for line in f.xreadlines():
        tabbed = line.replace('_', '\t')
        aln = tabbed.strip().split('\t')
        unique_reads[aln[0]] = [aln[1], str(chr[aln[2]]), aln[3], aln[5]]
    f.close()
    print str(len(unique_reads.keys())) + ' ' + read + ' reads align to S288c'
    # make SK1 chromosome notation compatible
    chr = {'chr01':1, 'chr02':2, 'chr03':3, 'chr04':4, 'chr05':5,
                   'chr06':6, 'chr07':7, 'chr08':8, 'chr09':9, 'chr10':10,
                   'chr11':11, 'chr12':12, 'chr13':13, 'chr14':14,
                   'chr15':15, 'chr16':16}
    # Read in SK1 alignments, line by line
    f = open(lib_name + '_SK1-' + fastq_num + '_' + read + '-data5.txt',
             'r')
    for line in f.xreadlines():
        tabbed = line.replace('_', '\t')
        aln = tabbed.strip().split('\t')
        if unique_reads.get(aln[0]):
            continue
        else:
            unique_reads[aln[0]] = [aln[1], str(chr[aln[2]]), aln[3], aln[5]]
    f.close()
    print str(len(unique_reads.keys())) + ' ' + read + \
               ' reads align to S288c or SK1'
    return unique_reads

def same_fragment(r1_dict, r2_dict, write_file=False):
    """extract reads that map to same restriction fragment in R1 and
    R2 reads. Option to write as tab-delimited file"""
    # make a dictionary with same fragment data where keys=Read_ID
    # items = Chr, Fragment, R1_Or, R1_Seq, R2_Or, R2_Seq 
    same_frag = {}
    for i in r1_dict:
        if r2_dict.get(i):
            if r1_dict[i][1] == r2_dict[i][1]:
                if r1_dict[i][2] == r2_dict[i][2]:
                    same_frag[i] = [r1_dict[i][1], r1_dict[i][2], r1_dict[i][0],
                                    r1_dict[i][3], r2_dict[i][0], r2_dict[i][3]]
    if write_file == True:
        print "need to add write file option... sorry!"
    return same_frag

def same_or(same_frag, lib_name, fastq_num):
    """extract reads that map to same restriction fragment in R1 and
    R2 reads in the same orientation. Writes output as tab-delimited file"""
    # extract entries from same_frag dictionary in same orientation in
    # R1 and R2
    same_frag_or = {}
    for i in same_frag:
        if same_frag[i][2] == same_frag[i][4]:
            same_frag_or[i] = '\t'.join(same_frag[i])
    lines = [i + '\t' + same_frag_or[i] for i in same_frag_or]
    f = open(lib_name + '_' + fastq_num + '-SFO.txt','w')
    f.write('\n'.join(i for i in lines) + '\n')
    f.close()
    return same_frag_or

def main(lib_name, fastq_num):
    # merge S288c and SK1 R1 read lists
    r1_dict = read_merge(lib_name, fastq_num, 'R1')
    # merge S288c and SK1 R2 read lists
    r2_dict = read_merge(lib_name, fastq_num, 'R2')
    # extract read pairs in same fragment
    same_frag = same_fragment(r1_dict, r2_dict)
    print str(len(same_frag.keys())) + " read pairs map to the same fragment"
    del r1_dict,r2_dict
    # extract read pairs in same orientation
    same_frag_or = same_or(same_frag, lib_name, fastq_num)
    print str(len(same_frag_or.keys())) + \
          " read pairs map to same fragment in the same orientation"


def main()
    pass

if __name__ == '__main__':
    pass