import argparse
import requests
import sys

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('-f', '--rettype',
                        choices=('fasta', 'docsum'),
                        default='fasta',
                        help='Return type')

    parser.add_argument('-i', '--accession-list',
                        help='List of accessions, one per line')

    parser.add_argument('-o', '--output',
                        help='Results file. Writes to stdout if not set')

    parser.add_argument('accessions',
                        nargs='*',
                        help='List of accession.version numbers')

    return parser.parse_args()

def gather_accessions(accession_list, accessions):
    '''Combines accessions from file accession_list (one acc per line)
    and the positional argument 'accessions'.

    Returns a list of all accessions.
    '''
    accs = []

    try:
        with open(accession_list, 'r') as f:
            for line in f:

                l = line.strip()

                accs.append(l)

    except FileNotFoundError:

        print('File {} does not exist!'.format(accession_list),
              file=sys.stderr)

        sys.exit(1)

    except TypeError:
        pass

    try:

        accs.extend(accessions)

    except TypeError:

        pass

    assert accs, 'No accessions have been provided'

    return accs

def fetch(rettype, accessions):
    '''Fetches requested accessions and returns them from NCBI'''

    acc_string = ' '.join(str(a) for a in accessions)

    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?\
            db=sequences&\
            id={accs}&\
            rettype={rettype}&\
            retmode=text'.replace(' ','').format(accs=acc_string,
                                                 rettype=rettype)

    res = requests.get(url)

    return res.text

def main():

    args = arguments()

    accs = gather_accessions(args.accession_list, args.accessions)

    fasta = fetch(args.rettype, accs)

    if args.output:

        with open(args.output, 'w') as o:
            o.write(fasta)

    else:

        print(fasta)

if __name__ == '__main__':
    main()
