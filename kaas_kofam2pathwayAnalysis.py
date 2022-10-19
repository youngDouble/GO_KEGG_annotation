#!/database/home/yangzhushuang/syssoft/python3/bin/python3
# -*- coding: utf-8 -*-

doc = """
####################################################################################################
Author:yangzhishuang  
E-mail:yangzs_chi@yeah.net
Edited by Yang Zhishuang at 2021/05/14
The latest version was edited at : 2021/05/14 By:yang zs 
####################################################################################################

This program is designed to perform kaas or KofamKOALA result. These files will be produced: *.kegg_anno.xls, *.kegg_pathway_stata.xls, *.K_codes_not_in_pathway_map.list.
Requirments:
    wget        windows: http://gnuwin32.sourceforge.net/packages/wget.htm
                unix: https://www.gnu.org/software/wget/

Usage:
python kaas_kofam2pathwayAnalysis.py  annot_kegg.txt [org]
    annot_kegg.txt      the first column should be the gene ID; Kd+ can be one or more in one row; the same geneID can be in one or more row.
                        Get from kaas or KofamKOALA result
    org      Organisms code of KEGG. eg. eco (Escherichia coli K-12 MG1655), The default value is "ko", Use "ko00001.keg" annotations not species-specific
             See: https://www.genome.jp/kegg/catalog/org_list.html

1. Use "ko00001.keg" annotations:
python kaas_kofam2pathwayAnalysis.py  annot_kegg.txt 

2. Use species-specific annotations:
python kaas_kofam2pathwayAnalysis.py  annot_kegg.txt eco


"""

import os
import sys
import re


def pathway_map(sp="ko"):
    """
    :param sp:The default is'ko', which means downloading 'ko00001.keg'
           https://www.genome.jp/kegg-bin/download_htext?htext=ko00001.keg&format=htext&filedir=
    :return: K_ko_map
    """
    #url = r"https://www.genome.jp/kegg-bin/download_htext?htext=asa00001.keg&format=htext&filedir="
    keg_file = sp +"00001.keg"
    cmd = r"wget 'http://www.kegg.jp/kegg-bin/download_htext?htext=" +sp + "00001.keg&format=htext&filedir=' -O " + keg_file
    if not os.path.exists(keg_file):
        try:
            res = os.system(cmd)
            # 使用system模块执行linux命令时，如果执行的命令没有返回值res的值是256
            # 如果执行的命令有返回值且成功执行，返回值是0
        except:
            print("Failed to run\n" + str(cmd) +"\nplease check the network")
            sys.exit()
    in_keg = open(keg_file, "r",encoding='gbk',errors='ignore').readlines()
    K_ko_map = {}
    for line in in_keg:
        if line.startswith("A"):
            # 'A09100 Metabolism'
            level_1 = re.match(r'^A(.+?)\s(.+)\n', line).group(2)
        elif line.startswith("B "):
            # B  09102 Energy metabolism
            level_2 = re.match(r'^B\s*(.+?)\s(.*)\n', line).group(2)
        elif line.startswith("C "):
            # 'C    00010 Glycolysis / Gluconeogenesis [PATH:asa00010]' or 'C    99980 Enzymes with EC numbers'
            pathway_info = re.match(r'^C\s*(\d+?)\s(.*)\n', line)
            pathway_id = "ko" + str(pathway_info.group(1))
            pathway_desc = str(pathway_info.group(2)).split(" [")[0]
            pathway_info_list = [pathway_desc, level_1, level_2]

        elif line.startswith("D ")   and in_keg[0] == '+D\tGENES\tKO\n':
            # 'D      ASA_1323 glk; glucokinase\tK00845 glk; glucokinase [EC:2.7.1.2]\n'
            K_info = re.match(r'^D\s*.*\t(K\d+?)\s(.*)\n', line)
            K_id = K_info.group(1)
            K_desc = K_info.group(2)
            K_info_list = [K_desc, pathway_id, pathway_desc, level_1, level_2]
            if K_id not in K_ko_map.keys():
                K_ko_map[K_id] = [K_info_list]
            else:
                l_tem = K_ko_map[K_id]
                if K_info_list not  in  l_tem:
                    l_tem.append(K_info_list)
                    K_ko_map[K_id] = l_tem

        elif line.startswith("D ")   and in_keg[0] == '+D\tKO\n':
            # 'D      K00844  HK; hexokinase [EC:2.7.1.1]'
            K_info = re.match(r'^D\s*(K\d+?)\s(.*)\n', line)
            K_id = K_info.group(1)
            K_desc = K_info.group(2)
            K_info_list = [K_desc, pathway_id, pathway_desc, level_1, level_2]
            if K_id not in K_ko_map.keys():
                K_ko_map[K_id] = [K_info_list]
            else:
                l_tem = K_ko_map[K_id]
                if K_info_list not in l_tem:
                    l_tem.append(K_info_list)
                    K_ko_map[K_id] = l_tem
    return (K_ko_map)

def ko_class_map():
    # https://www.genome.jp/kegg-bin/download_htext?htext=br08901.keg&format=htext&filedir=  htext
    # https://www.genome.jp/kegg-bin/download_htext?htext=br08901.keg&format=json&filedir=  josn
    # https://www.genome.jp/dbget-bin/get_linkdb?-t+orthology+path:ko00040
    cmd = r"wget 'https://www.genome.jp/kegg-bin/download_htext?htext=br08901.keg&format=htext&filedir=' -O br08901.keg"
    if not os.path.exists("br08901.keg"):
        try:
            res = os.system(cmd)
            # 使用system模块执行linux命令时，如果执行的命令没有返回值res的值是256
            # 如果执行的命令有返回值且成功执行，返回值是0
        except:
            print("Failed to run\n" + str(cmd) + "\nplease check the network")
            sys.exit()
    in_keg = open("br08901.keg", "r").readlines()
    pathway_map=open("kegg_pathway_map.xls","w+")
    ko_class = {}
    for line in in_keg:
        if line.startswith("A"):
            # 'A09100 Metabolism'
            level_1 = re.match(r'^A<b>(.*)</b>\n', line).group(1)
        elif line.startswith("B "):
            # B  09102 Energy metabolism
            level_2 = re.match(r'^B\s*(.*)\n', line).group(1)
        elif line.startswith("C "):
            # 'C    00010 Glycolysis / Gluconeogenesis [PATH:asa00010]' or 'C    99980 Enzymes with EC numbers'
            pathway_info = re.match(r'^C\s*(\d+?)\s(.*)\n', line)
            pathway_id = "ko" + str(pathway_info.group(1))
            pathway_desc = str(pathway_info.group(2)).split(" [")[0]
            line_out = pathway_id + "\t" + pathway_desc + "\t" + level_1 + "\t" +level_2 + "\t"
            pathway_map.write(line_out + "\n")
            pathway_info_list = [pathway_desc, level_1, level_2]
            if pathway_id not in ko_class.keys():
                ko_class[pathway_id] = pathway_info_list
    return(ko_class)


def K_list_Parser(K_list,sp="ko"):
    file_name = os.path.split(K_list)[1].rsplit(".",1)[0]
    out_kegg_anno = file_name + ".kegg_anno.xls"
    not_in_pathway_map = file_name + "K_codes_not_in_pathway_map.list"
    out_kegg_pathway = file_name + ".kegg_pathway_stata.xls"
    anno_f = open(out_kegg_anno, "w+")
    not_in_pathway_f = open(not_in_pathway_map, "w+")
    pathway_f = open(out_kegg_pathway, "w+")
    anno_f.write("gene_id\tK_id\tK_desc\tpathway_id\tpathway_desc\tlevel_1\tlevel_2\n")
    pathway_f.write("Pathway\tGenes annoted in term\tPathway ID\tLevel1\tLevel2\tKOs\tGenes\n")

    map_kegg = pathway_map(sp)
    ko_class = ko_class_map()
    infile_list = open(K_list, "r").readlines()
    infile_list = [ term.rstrip("\n").split("\t") for term in infile_list ]

    k_num_dict = {}
    line_out_tem = ""
    for line in infile_list:
        gene_id = line[0]
        if len(line) == 2:
            K_id = line[1]
            if K_id in map_kegg.keys():
                pathway_info = map_kegg[K_id]
                for K_info_list in pathway_info :
                    # K_desc, pathway_id, pathway_desc, level_1, level_2
                    string = "\t"
                    line_out =gene_id + "\t" + K_id + "\t" +string.join(K_info_list) + "\n"
                    if line_out != line_out_tem:
                        anno_f.write(line_out)
                        line_out_tem = line_out
            else:
                not_in_pathway_f.write(gene_id + "\t" + K_id + "\n")
                line_out = gene_id + "\t" * 6 + "\n"
                if line_out != line_out_tem:
                    anno_f.write(line_out)
                    line_out_tem = line_out
            # k_id 2 gene_id list
            if K_id not in k_num_dict.keys():
                k_num_dict[K_id] = gene_id
            else:
                k_num_dict[K_id] = k_num_dict[K_id] + ';' + gene_id
        else:
            anno_f.write(gene_id +"\t"*6 + "\n")

    ko_sample_dict = {}
    for K_id in list(k_num_dict.keys()):
        if K_id not in list(map_kegg.keys()):
            continue
        ko_num_sample = [ term[1] for term in  map_kegg[K_id]]
        for ko in ko_num_sample:
            if ko not in list(ko_sample_dict.keys()):
                ko_sample_dict[ko] = K_id
            else:
                ko_sample_dict[ko] = ko_sample_dict[ko] + ';' + K_id

    for ko in [item for item in list(ko_sample_dict.keys())  if  item  in list(ko_class.keys()) ] :
        pathway = ko_class[ko][0]
        level1 = ko_class[ko][1]
        level2 = ko_class[ko][2]
        k_num_list = ko_sample_dict[ko].split(';')
        gene_str = ''

        for k_num_sample in k_num_list:
            gene_str = gene_str + k_num_dict[k_num_sample] + ";"
            num_gene = gene_str.count(';')
        # pathway_f.write("Pathway\tGenes annoted in term\tPathway ID\tLevel1\tLevel2\tKOs\tGenes\n")
        pathway_f.write(pathway + '\t' + str(num_gene)+ '\t'  + ko +'\t' +level1 + '\t'+level2 +'\t'+ko_sample_dict[ko].rstrip(';')+ '\t' + gene_str.rstrip(';')+'\n')
    anno_f.close()
    pathway_f.close()



if len(sys.argv) < 2:  #直接执行本脚本给出帮助信息
    print(doc)
    sys.exit()
elif len(sys.argv) == 2:
        kaas_inflie = sys.argv[1]
        K_list_Parser(kaas_inflie)

elif len(sys.argv) == 3:
        kaas_inflie = sys.argv[1]
        sp = sys.argv[2]
        K_list_Parser(kaas_inflie,sp)
else:
    print(doc)
    sys.exit()

