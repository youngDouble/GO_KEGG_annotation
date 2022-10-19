#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author:yangzhishuang  
# E-mail:yangzs_chi@yeah.net
# Created on 2021/10/25 10:23
# The latest version was edited on : 2018/09/15  By:yangzhishuang


import sys
import os
import re
from argparse import ArgumentParser
from goatools import obo_parser, base
import pandas as pd


help_info = 'This script is convert the GO annotation with GO annotation with wego format, ' \
            'TSV format (interproscan.sh and eggNOG)  to GENE2TERM format, and '
def make_options():
    parser = ArgumentParser(description=help_info)
    parser.add_argument("-i", '--in', action='store', dest='inputfile',
                        help='path of input file',
                        default=False)

    avg = parser.parse_args()
    in_file = avg.inputfile
    # in_format = avg.format
    # return (in_file, in_format)
    return (in_file)



# in_file = r"G:\genebang\GO\TEST\test.go"
# main
def gene_go_query(goterm_list):
    '''
    :param goterm_list: eg. goterm_list = ['GO:0042026', 'GO:0006457', 'GO:0016887', 'GO:0005524']
    :return:
    '''
    all_goterm_list = goterm_list
    for goterm in goterm_list:
        goterm_obj = obofile.query_term(goterm)
        if goterm_obj != None:
            parents_list = list(goterm_obj.get_all_parents())
            all_goterm_list = list(set(all_goterm_list + parents_list))
        else:
            all_goterm_list.remove(goterm)
            continue
    # remove top level
    # GO:0003674 molecular_function
    # GO:0005575 cellular_component
    # GO:0008150 biological_process
    all_goterm_obj_list = [obofile.query_term(go_term) for go_term in all_goterm_list if go_term != "GO:0003674" and go_term != "GO:0005575" and go_term != "GO:0008150" ]
    all_goterm_detail_list = ['\t'.join([go_term.id, "level_" + str(go_term.level), go_term.namespace, go_term.name]) for go_term in all_goterm_obj_list if go_term is not None]
    return all_goterm_detail_list
    #

def  go_annotation_parser(in_file):
    '''
    :param in_anno:   Result of GO annotation with wego format,  TSV format of interproscan.sh and eggNOG.
    :return:
    '''
    in_anno_text = open(in_file, "r").readlines()
    in_anno_text = [term.rstrip("\n") for term in in_anno_text if term.startswith("#") == False]
    in_anno_GO = [[term.split("\t")[0], re.findall(r'GO:\d{7}', term)] for term in in_anno_text if len(re.findall(r'GO:\d{7}', term)) > 0]
    # in_anno_GO :
    # [['GENE1', ['GO:0003824']], ['GENE2', ['GO:0003824']], ['GENE3', ['GO:0005524', 'GO:0006457', 'GO:0016887'], ...]
    gene_GO_dict = {}
    for term in in_anno_GO:
        gene_id = term[0]
        GO_list = term[1]
        if gene_id  not in gene_GO_dict.keys():
            gene_GO_dict[gene_id] = GO_list
        else:
            gene_GO_dict[gene_id] = list(set(gene_GO_dict[gene_id] + GO_list))
    # gene_GO_dict :
    # {'GENE1': ['GO:0003824'], 'GENE2': ['GO:0008080'],'GENE3': ['GO:0042026', 'GO:0006457', 'GO:0005524'], ... }
    out_file.write('Gene_ID\tGO_term\tLevel\tFunction_class\tFunction\n')
    for gene_id in gene_GO_dict.keys():
        goterm_detail_list = gene_go_query(gene_GO_dict[gene_id])
        for GO_term in goterm_detail_list:
            out_file.write(gene_id + "\t" + GO_term + "\n")

    #To  DataFrame too slow
    # df_go = pd.DataFrame(columns=['Gene_ID', 'GO_term', 'Level', 'Function_class', 'Function'])
    # for gene_id in gene_GO_dict.keys():
    #     goterm_detail_list = gene_go_query(gene_GO_dict[gene_id])
    #     for GO_term in goterm_detail_list:
    #         df_go.loc[len(df_go)] = pd.Series({'Gene_ID': gene_id , 'GO_term': GO_term[0], 'Level': GO_term[1],
    #                                               'Function_class': GO_term[2], 'Function': GO_term[3]})
    # df_go.to_csv(out_file_path, sep='\t', index=False)



def r_plot(stat_file_name):
    content = """
library(ggplot2)
library(dplyr)
in_file <-  '""" + stat_file_name +"""'
out_file <- sub("_stats.xls","_stats_level2.pdf",in_file)

GO_table <- read.table(in_file,header = T,quote = '"', check.names = F,sep = "\\t")

colnames(GO_table) <-c("V1","V2","V3","V4")

GO_table$V5 <- paste(GO_table$V3,' (',GO_table$V1,')',sep = '')

# order tab
GO_table <- arrange(GO_table, desc(V2), desc(V4))
#GO_table <- arrange(GO_table, desc(V2), desc(V3))
level <- GO_table$V5
GO_table$V5<- factor(GO_table$V5, levels=level)


# color_list ==  pal_uchicago("default")(9)
#color_list <- c("#800000FF", "#767676FF", "#FFA319FF", "#8A9045FF", 
#                "#155F83FF", "#C16622FF", "#8F3931FF", "#58593FFF", 
#                "#350E20FF")
color_list <- c("#F0E685FF", "#466983FF", "#BA6338FF", "#5DB1DDFF")

pdf(file=out_file, bg="transparent",width = 12, height = 8) # save as pdf
par(mai=c(20,10,5,35),cex=6.7)

ggplot(GO_table, aes(x = V5 , y = V4 , fill = V2)) +
  # 条形图函数：stat表明取用样本点对应纵轴值
  # position_dodge(0.5) 表示同组内不同系列间错位比例 ：0.1表示90%重叠 
  geom_bar(stat = "identity" ,position = position_dodge(0.9),
           width =0.9,colour="black",size=0.2 )+
  # scale_color_manual(values = rep("black",26)) +
  scale_fill_manual(values= color_list)+
  # scale_fill_manual(values=brewer.pal(8, "Set2"))+
  geom_text(aes(label = V4 , y = V4 + 20 ), size=4,  position = position_dodge(0.5))+
  #去除绘图区和X轴之间的gap
  scale_y_continuous(expand = c(0, 0),limits = c(0,max(GO_table$V4)+80)) +
  labs(x = "GO Functional Classification", y = "No. of Genes")+
  coord_flip()+
  #theme_classic()
  theme_bw()+
  guides(color = guide_legend(ncol = 1),fill = guide_legend(reverse = F))+  #显示单列
  theme(legend.title=element_blank(), # 图例标题
        legend.background =element_blank(),# 图例背景
        legend.position = c("right"), #top,right,left, bottom
        legend.text = element_text( size = 9), #图例字体
        panel.grid = element_blank())

dev.off() 
"""
    with open("drawGO_level2.R", "w+") as drawGO_script:
        drawGO_script.write(content)


def go_anno_stat_and_drawing():
    '''
    Summarize and plot using GO annotations at the level 2
    '''
    df_go_anno = pd.read_table(out_file_path, header=0)
    df_go_anno_level2 = df_go_anno[df_go_anno[u'Level'] == "level_1"]
    df_level2_counts = df_go_anno_level2.groupby(["GO_term", "Function_class", "Function"], as_index=False)['Gene_ID'].count()
    df_level2_counts = df_level2_counts.rename(columns={'Gene_ID': 'No._of_Gene'})
    df_level2_counts.to_csv(stat_file_name, sep='\t', index=False)
    # write R script
    r_plot(stat_file_name)

# main over

if len(sys.argv) == 1:
    print(help_info)
    sys.exit()
else:
    in_file = make_options()

    # Download GO DAG file, go-basic.obo if not exist
    base.download_go_basic_obo("go-basic.obo")
    obofile = obo_parser.GODag('go-basic.obo', load_obsolete=False)

    path_list = os.path.split(in_file)
    out_file_name = path_list[-1].rsplit(".", 1)[0] + '_GO_anno.xls'
    out_file_path = os.path.join(path_list[0], out_file_name)
    stat_file_name = path_list[-1].rsplit(".", 1)[0] + '_GO_anno_stats.xls'
    stat_file_path = os.path.join(path_list[0], stat_file_name)

    out_file = open(out_file_path, "w+")
    # annotation
    go_annotation_parser(in_file)
    # stat and drawing
    go_anno_stat_and_drawing()
    os.system("Rscript ./drawGO_level2.R")
    os.remove("drawGO_level2.R")
    out_file.close()

