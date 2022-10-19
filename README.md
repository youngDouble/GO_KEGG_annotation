# GO_KEGG_annotation



## 01. kaas_kofam2pathwayAnalysis.py  

This script is designed to perform kaas or KofamKOALA result. These files will be produced: *.kegg_anno.xls, *.kegg_pathway_stata.xls, *.K_codes_not_in_pathway_map.list.

### **Requirments:**

- wget	windows: http://gnuwin32.sourceforge.net/packages/wget.htm, unix: https://www.gnu.org/software/wget/

### **Usage:**

```bash
python kaas_kofam2pathwayAnalysis.py  annot_kegg.txt [org]
    annot_kegg.txt      the first column should be the gene ID; Kd+ can be one or more in one row; the same geneID can be in one or more row.
                        Get from kaas or KofamKOALA result
    org      Organisms code of KEGG. eg. eco (Escherichia coli K-12 MG1655), The default value is "ko", Use "ko00001.keg" annotations not species-specific
             See: https://www.genome.jp/kegg/catalog/org_list.html

1. Use "ko00001.keg" annotations:
   python kaas_kofam2pathwayAnalysis.py  annot_kegg.txt 

2. Use species-specific annotations:
   python kaas_kofam2pathwayAnalysis.py  annot_kegg.txt eco
```

### Detail Information: 

[https://zhuanlan.zhihu.com/p/375740435?](https://zhuanlan.zhihu.com/p/375740435?)

### Cite:





## 02. GO_anno_from_tab.py

This script is convert the GO annotation with GO annotation with wego format, TSV format (interproscan.sh and eggnog-mapper)  to GENE2TERM format. (**Recommended interpro results**)

[**eggNOG**: http://eggnog5.embl.de/](http://eggnog5.embl.de/)

[**interpro**: https://www.ebi.ac.uk/interpro/search/sequence/](https://www.ebi.ac.uk/interpro/search/sequence/)

### **Requirments:**

```text
python3 (>=3.8.1)
    goatools （pip install goatools, Tested: goatools v1.1.6）
    pandas（pip install pandas, Tested: pandas v1.0.1  ）
R (>=3.6.2 )
    ggplot2 (>=3.3.5)
    dplyr (>=1.0.7)
```

### **Usage:**

```bash
python GO_anno_from_tab.py -i eggnog_result.emapper.annotations 
python GO_anno_from_tab.py -i go.wego
python GO_anno_from_tab.py -i interpro.out.txt
```

### Detail Information: 

[https://zhuanlan.zhihu.com/p/426424799?](https://zhuanlan.zhihu.com/p/426424799?)

###  Cite:

**goatools:** Klopfenstein D V, Zhang L, Pedersen B S, et al. GOATOOLS: A Python library for Gene Ontology analyses[J]. Scientific reports, 2018, 8(1): 1-17.

**ggplot2:** Valero-Mora P M. ggplot2: elegant graphics for data analysis[J]. Journal of Statistical Software, 2010, 35: 1-3.



