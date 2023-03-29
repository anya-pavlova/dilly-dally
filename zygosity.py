import argparse
from pathlib import Path
import os
import csv
import numpy as np
from collections import defaultdict
import pandas as pd
from IPython.display import display



def select_fields_tuple(d, keys):
    return tuple(d[key] for key in keys)


def select_fields(d, keys):
    return {key: d[key] for key in keys}


def is_nan(x):
    return x != x


def parse_vcf_info(info):
    result = {}
    for pair in info.split(";"):
        key_value = pair.split("=")
        if len(key_value) == 2:
            key, value = key_value
            if value == ".":
                value = np.nan
            result[key] = value
    return result

def read_vcf(vcf_file):
    with open(vcf_file) as f:
        while True:
            line = f.readline()
            if line.startswith("##"):
                continue
            elif line.startswith("#"):
                header = line.lstrip("#").split("\t")
                break
            else:
                raise Exception("Can't parse file")

        csv_reader = csv.reader(f, delimiter="\t")
        for row in csv_reader:
            row = dict(zip(header, row))
            row["BAM_INFO"] = row[header[-1]]
            row = {
                **select_fields(row, ("CHROM", "POS", "REF", "ALT", "QUAL", "BAM_INFO")),
                **parse_vcf_info(row["INFO"])
            }
            yield row


def count_homozygous_heterozygous_snv(path, sample_name, min_dp=None, min_qual=None):   
    #- 0/0 : the sample is homozygous reference
    #- 0/1 : the sample is heterozygous, carrying 1 copy of each of the REF and ALT alleles
    #- 1/1 : the sample is homozygous alternate
    chromosoms_counter = defaultdict(lambda: defaultdict(int))

    for row in read_vcf(path):
        if len(row["ALT"]) != 1:
            continue
        if 'DP' not in row:
            continue    

        if min_dp is not None and int(row["DP"]) < min_dp:
            continue
            
        if min_dp is not None and float(row["QUAL"]) < min_qual:
            continue
 
        if row["CHROM"] == "chrY" or row["CHROM"] == "chrX":
            continue
        
        chrom = row["CHROM"]
        zygocity = "homo" if "1/1" in row["BAM_INFO"] else "hetero"
        chromosoms_counter[chrom][zygocity] += 1
    
    flat_counter = [
        {"chr": key, "homo": value["homo"], "hetero": value["hetero"]}
        for key, value in chromosoms_counter.items()
    ]
    zygosity_table = pd.DataFrame(flat_counter)
    zygosity_table.set_index("chr")
    
    total_het_homo_zygosity = {'chr': 'total', 
                               'homo': zygosity_table['homo'].sum(), 
                               'hetero': zygosity_table['hetero'].sum()}
    
    zygosity_ratio = {'chr': 'het/homo ratio',
                      'homo': zygosity_table['hetero'].sum() / zygosity_table['homo'].sum(),
                      'hetero': ' '}
    
    zygosity_table = zygosity_table.append(total_het_homo_zygosity, ignore_index=True)
    zygosity_table = zygosity_table.append(zygosity_ratio, ignore_index=True)

    zygosity_table['homo'] = (zygosity_table['homo'] * 10).astype(np.int64) / 10
    
    zygosity_table.to_csv(path.parent / (sample_name + '_zygosity_level.csv'), index=False)   
    
    return zygosity_table


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate zygosity from vcf')
    parser.add_argument('vcf_path', metavar='N', help='Path to sample directory')
    args = parser.parse_args()

    vcf_path = Path(args.vcf_path)
    sample_name = vcf_path.name
    
    for file in vcf_path.iterdir():
        if file.name.endswith("norm_variants_only_by_target.vcf") and file.name.startswith('merged_'):
            zygosity = count_homozygous_heterozygous_snv(file, sample_name, min_dp = 13, min_qual = 0)

    



    


