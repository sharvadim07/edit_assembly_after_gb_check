import logging
import argparse
from Bio import SeqIO
from Bio.Seq import Seq


parser = argparse.ArgumentParser(description='Program for editing assembly after genbank contamination check.')
parser.add_argument('-a','--assembly', type=str,                   
                    help='Input assembly file.', required=True)
parser.add_argument('-c','--check_results', type=str,                   
                    help='File with genbank checking report.', required=True)

args = parser.parse_args()

def read_genbank_report(report_file_name):
    with open(report_file_name, 'r') as report_file:
        exclude_section = False
        trim_section = False
        duplicate_section = False
        exlude_seq_list = []
        trim_section_dict = {}
        duplicate_list = []
        for report_line in report_file:
            if report_line == '\n':
                exclude_section = False
                trim_section = False
                duplicate_section = False                        
            elif report_line.strip() == 'Exclude:':
                exclude_section = True
                trim_section = False
                duplicate_section = False
            elif exclude_section and 'Sequence name, length, apparent source' not in report_line:
                exlude_seq_list.append(report_line.strip().split()[0])
            elif report_line.strip() == 'Trim:':
                exclude_section = False
                trim_section = True
                duplicate_section = False
            elif trim_section and 'Sequence name, length, span(s), apparent source' not in report_line:
                report_line_splitted = report_line.strip().split()
                seq_name = report_line_splitted[0]
                span_list = report_line_splitted[2].split(',')
                trim_section_dict[seq_name] = span_list
            elif report_line.strip() == 'Duplicated:':
                exclude_section = False
                trim_section = False
                duplicate_section = True
            elif duplicate_section and 'Sequence names, length' not in report_line:
                duplicate_list.append(report_line.strip().split()[1].split('|')[1])
    return exlude_seq_list, trim_section_dict, duplicate_list


def split_by_list_seq(in_seq, separator_seq_list):    
    import re
    return re.split(r'(?:' + '|'.join(separator_seq_list) + r')', in_seq)

def edit_assembly(fasta_file_name, exlude_seq_list, trim_section_dict, duplicate_list):
    assembly_dict = SeqIO.to_dict(SeqIO.parse(fasta_file_name, 'fasta'))
    for seq in exlude_seq_list + duplicate_list:
        if seq in assembly_dict:
            assembly_dict.pop(seq)
    cont_seq_dict = {}
    for seq in trim_section_dict:
        if seq in assembly_dict:
            for span in trim_section_dict[seq]:
                start_pos = int(span.split('..')[0]) - 1
                end_pos = int(span.split('..')[1])
                cont_seq_dict[seq + '_' + span] = assembly_dict[seq].seq._data[start_pos:end_pos]
            #TODO: check this
            new_seq_list = split_by_list_seq(assembly_dict[seq].seq._data, [cont_seq_dict[cont_seq] for cont_seq in cont_seq_dict])
            for i, new_seq in enumerate(new_seq_list):
                if new_seq != '':
                    name = seq + '_trimmed_' + str(i)
                    assembly_dict[name] = SeqIO.SeqRecord(Seq(new_seq), name, name, name)
            assembly_dict.pop(seq)
    with open('contamination.fasta', 'w') as contamination_file:
        for cont_seq in cont_seq_dict:
            contamination_file.write('>' + cont_seq + '\n')
            contamination_file.write(cont_seq_dict[cont_seq] + '\n')


    import os
    with open('corr_' + os.path.basename(fasta_file_name), 'w') as handle:
        SeqIO.write(assembly_dict.values(), handle, 'fasta')


exlude_seq_list, trim_section_dict, duplicate_list = read_genbank_report(args.check_results)
edit_assembly(args.assembly, exlude_seq_list, trim_section_dict, duplicate_list)

