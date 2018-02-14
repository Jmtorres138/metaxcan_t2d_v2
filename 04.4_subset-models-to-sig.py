#!/usr/bin/python -O
# Jason Matthew Torres
'''
Subset full model snps file to only those tissue models identified as significant
Usage:
'''
# libraries
import sys,os
import gzip,time

cur_dir = "/well/mccarthy/users/jason/projects/metaxcan_t2d_v2/"
df_dir = cur_dir + "data_frames/"
model_file = df_dir + "model-snps-v6p_full.txt.gz"

gws_sig_file = df_dir + "Scott_T2D_gws-significant.full.txt"
lws_sig_file = df_dir + "Scott_T2D_lws-significant.full.txt"

def build_sig_dic(sig_file):
    fin = open(sig_file,'r')
    dic = {}
    header = fin.readline()
    for line in fin:
        l = line.strip().split()
        ensid, best_tissue, tissue = l[0], l[5], l[11]
        if best_tissue=="NA":
            tiss = tissue
            if tiss=="Muscle-skeletal":
                tiss="Muscle_Skeletal"
            if tiss=="Adipose-subcutaneous":
                tiss="Adipose_Subcutaneous"
        else:
            tiss = best_tissue
        key = ensid+":"+tiss
        dic[key] = key
        #try:
        #    dic[ensid].append(tiss)
        #except:
        #    dic[ensid] = [tiss]
    fin.close()
    return dic


def subset_and_write_modsnps(gws_dic,lws_dic,gws_out_name,lws_out_name):
    fin = gzip.open(model_file,'rb')
    fout1 = gzip.open(gws_out_name,'wb')
    fout2 = gzip.open(lws_out_name,'wb')
    head_list = fin.readline().strip().split()
    fout1.write("\t".join(head_list)+"\n")
    fout2.write("\t".join(head_list)+"\n")
    count=0
    for line in fin:
        count+=1
        sys.stdout.write("\r%d" % count)
        sys.stdout.flush()
        l = line.strip().split()
        ensid,rsid,tiss,chrom,pos = l[0],l[1],l[2],l[3],l[4]
        key = ensid+":"+tiss
        try:
            check = len(gws_dic[key])
            fout1.write("\t".join(l)+"\n")
            check = len(lws_dic[key])
            fout2.write("\t".join(l)+"\n")
        except:
            try:
                check = len(gws_dic[key])
                fout1.write("\t".join(l)+"\n")
            except:
                try:
                    check = len(lws_dic[key])
                    fout2.write("\t".join(l)+"\n")
                except:
                    pass
    fin.close()
    fout1.close()
    fout2.close()




def main():
    gws_dic =  build_sig_dic(gws_sig_file)
    lws_dic =  build_sig_dic(lws_sig_file)
    subset_and_write_modsnps(gws_dic,lws_dic,df_dir+"model-snps-v6p_full_gws.txt.gz",df_dir+"model-snps-v6p_full_lws.txt.gz")

if (__name__=="__main__"): main()
