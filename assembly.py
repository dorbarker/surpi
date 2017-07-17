import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory
from Bio import SeqIO
from utilities import concatenate, fastq_to_fasta

def fastq_seq_length(inputfile: Path) -> int:
    '''Returns the length of the first sequence in a FASTQ file'''

    with inputfile.open('r') as fastq:
        for rec in SeqIO.parse(fastq, 'fastq'):
            length = len(rec.seq)
            break
        return length

def sequniq(matched: Path, uniqvir: Path):
    '''Extracts the unique sequences from `matched`'''

    sequniq_cmd = ('gt', 'sequniq', '-seqit', '-force',
                   '-o', str(uniqvir), str(matched))

    subprocess.check_call(sequniq_cmd)

def abyss(unmatched_add_vir: Path, length: int, cores: int, kmer: int,
          ignore_barcodes: bool):
    '''Executes abyss_minimus.sh'''

    contig_cutoff = int(1.75 * length)

    abyss_cmd = map(str, ('abyss_minimus.sh', unmatched_add_vir, length,
                          contig_cutoff, cores, kmer, ignore_barcodes))

    subprocess.check_call(abyss_cmd)

def assemble(matched_vir_fastq: Path, uniqunmatched: Path, workdir: Path,
             tempdir: Path, sample_fastq: Path,
             kmer: int, ignore_barcodes: bool, cores: int) -> Path:

    # matched_vir_fastq is temporary
    # uniqunmatched comes from extract_to_fast
    # matched_vir_uniq is created by sequniq, only used here,
        # then sent to workdir
    # addvir is created by sequniq
    # sample_fastq_length comes externally?

    sample_fq_length = fastq_seq_length(sample_fastq)
    contig_cutoff = int(1.75 * sample_fq_length)

    addvir = (workdir / matched_vir_fastq.stem).with_suffix('.addVir.fasta')

    # TODO: change horrible name
    output_name = 'all.{}.unitigs.cut{}.{}-mini.fasta'.format(addvir.name,
                                                              sample_fq_length,
                                                              contig_cutoff)

    output = addvir.with_name(output_name)

    with TemporaryDirectory(dir=str(tempdir)) as ephemeral:
        throw_away_base = Path(ephemeral) / matched_vir_fastq.stem

        matched_vir_fasta = throw_away_base.with_suffix('.fasta')
        matched_vir_uniq = throw_away_base.with_suffix('.uniq.fasta')

        fastq_to_fasta(matched_vir_fastq, matched_vir_fasta)

        sequniq(matched_vir_fasta, matched_vir_uniq)

        concatenate(matched_vir_uniq, uniqunmatched, output=addvir)

    abyss(addvir, sample_fq_length, cores, kmer, ignore_barcodes)

    return output
