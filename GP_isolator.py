import time
from time import strftime

#opens the txt file containing nucleotide bases
#collects contents of file 
with open('raw data\ebola_nucleotide_bases.txt') as raw_fasta_file:
    fasta_content = [i.replace('\n','') for i in raw_fasta_file]

#collects acession codes - line that containts >
acession_codes = [line_no for line_no in fasta_content if '>' in line_no]

#collects the index acession codes in fasta_content
acession_code_index = [data[0] for data in enumerate(fasta_content) if '>' in data[1]]

#generates the range values (relative to fasta_content) of each acession code's nucleotide data
acession_code_ranges = []
for count, indx in enumerate(acession_code_index):    
    try:
        acession_code_ranges.append([indx+1, acession_code_index[count+1]])
    except IndexError:
        acession_code_ranges.append([indx+1, len(fasta_content)])
        
#collects nucleotide data of each acession code based on the generated range values
acession_code_bp = []
for line_range in acession_code_ranges:
    indiv_bp = []
    data_lines = [lines for lines in fasta_content[line_range[0]:line_range[1]]]    
    
    for line in data_lines:
        for bp in line:
            indiv_bp.append(bp)
    
    acession_code_bp.append(indiv_bp)

#opens txt file containing the locations of glycoprotein gene in each acession code
#generates range values glycoprotein nucleotide bases relative to each acession code's nucleotide data
gp_ranges = []
with open('raw data\glycoprotein_location.txt') as gp_loc_file:
    gp_locs = [line.replace('\n','') for line in gp_loc_file]
    
    for gene_range in gp_locs:
    
        start = int(gene_range.split('..')[0])
        end = int(gene_range.split('..')[1])
        
        gp_ranges.append([start,end])

#collects the actual nucleotide bases of the glycoprotein gene on each acession code, based on the generated glycoprotein range values
acession_code_gp_data = []
for indx, range_val in enumerate(gp_ranges):
    gp_data = ''
    for bp in acession_code_bp[indx][range_val[0]-1:range_val[1]]:
        gp_data += bp
        if len(gp_data) % 70 == 0:
            gp_data += '\n'
        else:
            continue
        
    acession_code_gp_data.append(gp_data)


#time_stamp = strftime("%b %d %Y") could be used to document time of file generation
#generates a file containing the acession code, and name of strain/species; and its isolated glycoprotein nucleotide bases
with open(f'results\ebola_isolated_gp.txt', 'w') as raw_fasta_file:
    for count, code in enumerate(acession_codes):    
        raw_fasta_file.write(f'>{code}\n')
        for bp in acession_code_gp_data[count]:
            raw_fasta_file.write(bp)
        raw_fasta_file.write('\n\n')

print('=file written=')