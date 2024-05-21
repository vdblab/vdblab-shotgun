import sys
import subprocess
import pandas as pd

main_dir = '/data/peledj/baichom1/Projects/Peled_Analysis/2024/01_Enterococcus_Nutrition_MTX/cazyme_pipeline/'
num_reads_file = ''.join([main_dir, 'num_reads_per_exp.txt'])
dbcan_dir = ''.join([main_dir, 'dbcan/'])
align_dir = ''.join([main_dir, 'bwa_align/'])
save_dir = ''.join([main_dir, 'detailed_dbcan/'])

def write_dbcan_info_file(aligned_file, substrate_file, cgc_file, overview_file, save_file, num_reads):

    overview_df = pd.read_csv(overview_file, sep = "\t")
    counts_df = pd.read_csv(aligned_file, sep = '\t', names=['bam_name', 'zero', 'length', 'counts']).drop('zero', axis = 'columns')
    
    # from the substrate file we will take the predicted substrate and the "k_n" name they are assigned to:
    s_k_n = []; s_sub = []
    with open(substrate_file, 'r') as sub_f:
        line = sub_f.readline()
        while line:
            if line.startswith('k'):
                tab_parts = line.split('\t')
                s_k_n.append(tab_parts[0].split('|')[0])
                s_sub.append(tab_parts[2])
            line = sub_f.readline()
    sub_df = pd.DataFrame({
        'k_name' : s_k_n,
        'substrate' : s_sub,
    })

    # the cgc file acts as a "key file" connecting the "k_n"s used in the substrate file, and the bam_gene names;
    bam_n = []; k_n = []; type_n = []
    with open(cgc_file, 'r') as cgc_f:
        line = cgc_f.readline()
        while line:
            if not line.startswith('+'):
                sub_lines = line.split("\t")
                bam_n.append(sub_lines[8])
                k_n.append(sub_lines[5])
                type_n.append(sub_lines[1])
            line = cgc_f.readline()
    cgc_faa_df = pd.DataFrame({
        'bam_name' : bam_n,
        'k_name' : k_n,
        'caz_type' : type_n,
    })


    caz_df = overview_df.merge(counts_df, how = 'left', left_on = "Gene ID", right_on = 'bam_name').drop('bam_name', axis= 'columns')
    caz_df = caz_df.merge(cgc_faa_df, how = 'left', left_on = "Gene ID", right_on = "bam_name")
    caz_df = caz_df.merge(sub_df, how = 'left', on = 'k_name')
    caz_df['RPM'] = ((10**6)*caz_df['counts'])/num_reads
    caz_df.to_csv(save_file, sep = '\t')

def count_num_reads_compressed_file(file_name):
    ''' 
    Small helper function which will take a compressed fastq file and return the number of reads.
        -file_name = str (or path) the location of the file to count number of reads for. 
    '''
    command = f"echo $(($(zcat {file_name} | wc -l)/4))"
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE)
    return float(result.stdout)

if __name__ == "__main__":
    '''
    Arguments required (given in the order expected:):
        - aligned_file = the depth file generated from aligning the host-depleted 
        reads to the metaerg annotated genes. 
        - substrate_file = the run_dbcan substrate file. 
        - cgc_file = the run_dbcan cgc.out file. 
        - overview_file = the run_dbcan overview file. 
        - save_file = name of the file to write the RPM outputs to.
        - r1 = the location of the r1 files for calculating the number of reads in the host-depleted file. 
    '''
    if "snakemake" not in globals():
        # assume taking in variables from the main argv:
        aligned_file = sys.argv[1]
        substrate_file = sys.argv[2]
        cgc_file = sys.argv[3]
        overview_file = sys.argv[4]
        save_file = sys.argv[5]
        r1 = sys.argv[6]
    else:
        aligned_file = snakemake.input.coverage
        substrate_file = snakemake.input.substrate
        cgc_file = snakemake.input.cgc
        overview_file = snakemake.input.overview
        save_file = snakemake.output.rpm_file
        r1 = snakemake.input.r1
    num_reads = count_num_reads_compressed_file(str(r1))
    print(f"Num reads: {num_reads}")
    write_dbcan_info_file(aligned_file, substrate_file, cgc_file, overview_file, save_file, num_reads)