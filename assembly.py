import shutil
import subprocess
import tempfile
from typing import List, Dict
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
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

def abyss(unmatched_addvir: Path, length: int, cores: int, kmer: int,
          ignore_barcodes: bool, workdir: Path, tempdir: Path):
    '''Executes abyss_minimus.sh'''

    contig_cutoff = int(1.75 * length)

    outputname = '{}.unitigs.cut{}.{}-mini.fasta'.format(unmatched_addvir.name,
                                                         length, contig_cutoff)
    out = workdir / outputname

    reads_by_barcodes = demultiplex(unmatched_addvir, ignore_barcodes)

    to_cat = []

    with tempfile.TemporaryDirectory(dir=str(tempdir)) as ephemeral:
        for barcode, reads in reads_by_barcodes.items():

            unitigs = run_abyss(reads, kmer, barcode, Path(ephemeral), cores)
            adjust_fasta(unitigs, kmer, barcode, length)

            mini = minimo(unitigs)

            adjust_fasta(mini, kmer, barcode, contig_cutoff)

            new_name = workdir / mini.name
            shutil.move(str(mini), str(new_name))

            to_cat.append(new_name)

    concatenate(*to_cat, output=out)

    return out

def assemble(matched_vir_fastq: Path, uniqunmatched: Path, workdir: Path,
             tempdir: Path, sample_fastq: Path,
             kmer: int, ignore_barcodes: bool, cores: int) -> Path:

    sample_fq_length = fastq_seq_length(sample_fastq)

    addvir = (workdir / matched_vir_fastq.stem).with_suffix('.addVir.fasta')


    with tempfile.TemporaryDirectory(dir=str(tempdir)) as ephemeral:
        throw_away_base = Path(ephemeral) / matched_vir_fastq.stem

        matched_vir_fasta = throw_away_base.with_suffix('.fasta')
        matched_vir_uniq = throw_away_base.with_suffix('.uniq.fasta')

        fastq_to_fasta(matched_vir_fastq, matched_vir_fasta)

        sequniq(matched_vir_fasta, matched_vir_uniq)

        concatenate(matched_vir_uniq, uniqunmatched, output=addvir)

    output = abyss(addvir, sample_fq_length, cores, kmer,
                   ignore_barcodes, workdir, tempdir)

    return output

def demultiplex(fastafile: Path, ignore_barcodes: bool):

    tr = str.maketrans('#/', '\t\t')

    reads_by_barcodes = {}  # type: Dict[str, SeqRecord]
    with fastafile.open('r') as fasta:
        for rec in SeqIO.parse(fasta, 'fasta'):

            if ignore_barcodes:
                try:
                    reads_by_barcodes['all_reads'].append(rec)
                except KeyError:
                    reads_by_barcodes['all_reads'] = [rec]

            else:

                _, barcode, *_ = rec.id.translate(tr).split('\t')

                try:
                    reads_by_barcodes[barcode].append(rec)

                except KeyError:
                    reads_by_barcodes[barcode] = [rec]

    return reads_by_barcodes

def run_abyss(fasta_reads: List[SeqRecord], kmer: int, name: str, workdir: Path, cores: int) -> Path:

    fastafile = (workdir / name).with_suffix('.fasta')
    unitigs = workdir / (name + '-unitigs.fa')

    with fastafile.open('w') as fasta:
        SeqIO.write(fasta_reads, fasta, 'fasta')

    cmd = ('abyss-pe', 'k={}'.format(kmer), 'name={}'.format(workdir / name),
           'se={}'.format(fastafile), 'np={}'.format(cores))

    subprocess.check_output(cmd, stderr=subprocess.STDOUT)

    return unitigs

def minimo(abyss_fasta: Path) -> Path:

    cmd = ('Minimo', str(abyss_fasta), '-D', 'FASTA_EXP=1')
    subprocess.check_output(cmd, stderr=subprocess.STDOUT)

    return abyss_fasta.parent / (abyss_fasta.stem + '-contigs.fa')

def adjust_fasta(fasta: Path, kmer: int, barcode: str, length_cutoff: int):

    outsequences = []

    with fasta.open('r') as infasta:

        for rec in SeqIO.parse(infasta, 'fasta'):

            new_description = '{}_kmer{}_barcode{}'.format(rec.description,
                                                           kmer, barcode)

            rec.description = ''

            if barcode not in rec.id:
                rec.id = '_'.join((rec.id, new_description.replace(' ', '_')))

            if len(rec.seq) >= length_cutoff:
                outsequences.append('>{}\n{}'.format(rec.id, rec.seq))


    # oneline sequences, since some downstream work requires this
    fasta.write_text('\n'.join(outsequences))
