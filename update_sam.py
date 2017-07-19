from pathlib import Path

def is_sam_header(line: str) -> bool:
    '''If the line begins with "@" return True, else False'''
    return line.startswith('@')

def compare_sam(infile: Path, outfile: Path):

    new_hits, existing_hits = 0, 0
    total_edit_dist_new, total_edit_dist_existing = 0, 0
    replacements = 0
    sum_data = 0

    with infile.open('r') as annot, outfile.open('w') as output:

        for line in annot:

            data = line.split()

            if is_sam_header(line):

                output.write(line)
                continue

            first_entry = data[0].split('|')

            edit_distance = int(data[12].split(':')[2])

            try:
                edit_distance_prev = int(first_entry[1])
            except IndexError:
                edit_distance_prev = 0

            if edit_distance >= 0:

                new_hits += 1
                total_edit_dist_new += edit_distance

                if len(first_entry) > 1:  # hit

                    existing_hits += 1
                    total_edit_dist_existing += edit_distance_prev

                    if edit_distance <= edit_distance_prev:

                        replacements += 1
                        sum_data += edit_distance

                        acc = data[2]

                        original = '|{}|{}'.format(*first_entry[:2])
                        rep = '|{}|{}'.format(edit_distance, acc)

                        output_line = line.replace(original, rep, 1)

                        output.write(output_line)

                    else:

                        sum_data += edit_distance_prev
                        output_line = line

                else:

                    existing_hits += 1
                    sum_data += edit_distance
                    acc = data[2]

                    rep = '{}|{}|{}'.format(data[0], edit_distance, acc)

                    output_line = line.replace(data[0], rep, 1)

                    replacements += 1

            else:

                output_line = line

                if len(first_entry) > 1:

                    existing_hits += 1

                    sum_data += edit_distance_prev

                    total_edit_dist_existing += edit_distance_prev

            output.write(output_line)

def update_sam(infile: Path, outfile: Path):

    with infile.open('r') as annot, outfile.open('w') as output:

        for line in annot:

            if is_sam_header(line):

                output.write(line)
                continue

            data = line.split()

            header = data[0].split('|')

            if len(header) > 1:

                fst, dvalue, acc, *_ = header

                edit_distance = int(data[12].split(':')[2])

                output_line = line.replace(data[0], fst)

                if edit_distance >= 0:

                    rep = 'NM:i:{}'.format(dvalue)
                    output_line = output_line.replace(data[13], rep)

                else:
                    rep = '{}\tNM:i:{}'.format(data[12], dvalue)
                    output_line = output_line.replace(data[12], rep)

                output.write(output_line)

            else:

                output.write(line)
