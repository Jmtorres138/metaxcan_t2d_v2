#!/usr/bin/python -O
# Jason Matthew Torres
'''
Build key with CHR and POS information for rsids where available
Usage: python JTbuildSNPkey-fdr05.py
'''
# libraries
import sys,os
import gzip,time

cur_dir = "/well/mccarthy/users/jason/projects/metaxcan_t2d_v2/model_snps/GTEx-V6p-HapMap-2016-09-08/"
model_file = cur_dir + "model-snps-v6p.txt"
ref_dir = "/well/got2d/jason/reference/dbSNP147-hg19-allSNPs/"
ref_file = ref_dir + "dbSNP147.hg19.allSNPs.bed.gz"


def snp_list_dic():
    print ("Building dictionary of gwas RSIDS...")
    fin = open(model_file,'r')
    fin.readline()
    snplist = []
    dic = {}
    count=0
    for line in fin:
        count+=1
        sys.stdout.write("\r%d" % count)
        sys.stdout.flush()
        l = line.strip().split()
        snpid = l[1]
        dic[snpid] = snpid
        snplist.append(snpid)
    fin.close()
    return(snplist,dic)


def write_keyfile(snplist, dic):
    out_name = cur_dir + "snp_keyfile.txt.gz"
    fin = gzip.open(ref_file,'rb')
    fin.readline()
    ref_dic = {}
    print "\nBuilding referencence dictionary..."
    count=0
    for line in fin:
        count+=1
        sys.stdout.write("\r%d" % count)
        sys.stdout.flush()
        l = line.strip().split()
        rsid, chrom, pos = l[3], l[0], l[2]
        try:
            len(dic[rsid])
            ref_dic[rsid] = [rsid,chrom,pos]
        except:
            pass
    fin.close()
    print ("\nWriting key file..")
    fout = gzip.open(out_name,'wb')
    head_list = ["RSID","CHR","POS"]
    fout.write("\t".join(head_list)+"\n")
    count=0
    for snp in snplist:
        count+=1
        sys.stdout.write("\r%d" % count)
        sys.stdout.flush()
        try:
            write_list = ref_dic[snp]
        except:
            try:
                temp_list = snp.split(":")
                chrom,pos = temp_list[0], temp_list[1]
                write_list = [snp,chrom,pos]
            except:
                write_list = [snp, "NA","NA"]
        fout.write("\t".join(write_list)+"\n")
    fout.close()
    print("\nProcess Complete")


def main():
    snplist, dic = snp_list_dic()
    write_keyfile(snplist,dic)


if (__name__=="__main__"): main()
