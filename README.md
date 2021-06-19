

# tslysis
***

tslysis is a R package that includes scripts and functions to use extracellular rRNA as the proxy to estimate the taxon-specific microbial lysis, lysis rate, mortality rate and etc in water. Previously, Zhong et al. (2021) reported the presence of millions copies of extracellular rRNA in one-millimetre seawater. These extracellular free rRNAs, probably in the form of ribonucleoprotein-bound rRNAs such as ribosomes, are shown to be the result of cell lysis (e.g. viral lysis) and can be used to estimate the taxon-specific lysis in the sea, which provides an important tool in our quest to understand the distribution and abundance of microbes in nature.


&nbsp;
&nbsp;


![](https://github.com/kevinzhongxu/CasOligo/tree/master/inst/extdata/oligo_distribution_among_taxa.png)
**Fig.1** Distribution of the extracellular rRNA in costal seawater at Strait of Geogia nearby Quadra, BC, Canada. 


&nbsp;
&nbsp;


# Features
***
  * The lysis() function allows calculate taxon-specific lysis using the ratio of the relative abundance of extracellular rRNA to cellular rRNA.
  
  * The lysis_rate() function allows estimate taxon-specific lysis rate. It is calculated using the extracellular rRNA turnover rate by multiplying the ratio of the concentration of extracellular rRNA to cellular rRNA.
  
  * ......
  
&nbsp;
&nbsp;


# Installation
***

To install the latest version from GitHub, simply run the following from an R console:

```r
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("kevinzhongxu/tslysis")
```

&nbsp;
&nbsp;
&nbsp;

# Dependancy
***
This package depends on the pre-installation of following R package: 

  * VennDiagram
  * reshape2
  
&nbsp;
&nbsp;


# Citation
***
&nbsp;

If you use tslysis in a publication, please cite our article in [here](https://www.biorxiv.org/content/10.1101/2020.06.02.130807v1):

Zhong KX, Wirth JF, Chan AM & Suttle CA. Extracellular ribosomal RNA provides a window into taxon-specific microbial lysis in the sea. xxx xx(xx): xxxx (2021). https://www.biorxiv.org/content/10.1101/2020.06.02.130807v1

&nbsp;
&nbsp;

# Get start
***
&nbsp;

## Example 1: estimate taxon-specific lysis 
This is an example to estimate taxon-specific lysis


```r
# If you want to estimate the taxon-specific lysis for prokaryotes, then run
lysis(otu_table="Path/To/Your/OTU_TABLE_FILE_16S.txt")


# If you want to test the funtion/script, then you can run lysis() function using one otu table example from the package:
otu_table <- read.table(as.character(system.file("extdata", "quadra_16S.txt", package = "tslysis")), h=T, sep="\t", quote=NULL, comment='', fill=T, stringsAsFactors = F)

lysis(otu_table=otu_table)

```

The input file is a tabular-separated otu table (.txt) contains one column named with "OTU_ID", another column named with "taxonomy", and the remaining columns named by "sample.id"_"one of three rDNA/rRNA fraction".  For example, if a water sample with "sample.id"="May_30m", then in this otu_table the column names for each rDNA/rRNA fraction of this sample is suggested to be: May_30m_Cellular.rDNA, May_30m_Cellular.rRNA and May_30m_Extracellular.rRNA.


```r
# An example of out_table input for lysis()
otu_table <- read.table(as.character(system.file("extdata", "example_out_table.txt", package = "tslysis")), h=T, sep="\t", quote=NULL, comment='', fill=T, stringsAsFactors = F)

print(otu_table)

```

lysis.R function will firstly calculate the relative abundance of each fraction of the samples. Then the lysis will be calculated by the ratio of the relative abundance of extracellular.rRNA to cellular.rRNA. Thus, it's important to be noted that the taxa within the otu_table should be pretreated to include only the microbial group the you intend to study. For example, if you aim for calculate the taxon-specific lysis for prokaryotes, then the "Eukaryota", "Chloroplast" and "Mitochondria" should be removed from the otu table. Another example, if you aim for calculate the taxon-specific lysis for microeukaryotes then the "Metazoa", "Bacteria" and "Archaea" should be removed from the otu table

&nbsp;


## Example 2: estimate taxon-specific lysis rate
This is an example to estimate taxon-specific lysis rate


```r
# If you want to estimate the taxon-specific lysis rate for prokaryotes, then run
lysis_rate(otu_table="Path/To/Your/OTU_TABLE_FILE_16S.txt", metadata="Path/To/Your/METADATA.txt")


# If you want to test the funtion/script, then you can run lysis_rate() function using the otu_table and metadata example from the package:
otu_table <- read.table(as.character(system.file("extdata", "quadra_16S.txt", package = "tslysis")), h=T, sep="\t", quote=NULL, comment='', fill=T, stringsAsFactors = F)

metadata <- read.table(as.character(system.file("extdata", "metadata.txt", package = "tslysis")), h=T, sep="\t", quote=NULL, comment='', fill=T, stringsAsFactors = F)

lysis_rate(otu_table=otu_table, metadata=metadata)

```

Compared to lysis(), lysis_rate() requires another input file, which is a tabular-separated metadata file (.txt) contains at least columns "sample.id", "Abundance_Total_Cellular.rRNA", "Abundance_Total_Extracellular.rRNA" and "turnover_rate". 

For example, 
```r
# An example of metadata input for lysis_rate()
metadata <- read.table(as.character(system.file("extdata", "metadata.txt", package = "tslysis")), h=T, sep="\t", quote=NULL, comment='', fill=T, stringsAsFactors = F)

print(metadata)

```

"sample.id" = "May_30m", which is the name of water sample.
"Abundance_Total_Cellular.rRNA" = Abundance of total cellular rRNA measured from the water sample.
"Abundance_Total_Extracellular.rRNA" = Abundance of total extracellular rRNA measured from the water sample.
"turnover_rate" = The estimated turnover rate of extracellular rRNA in the water sample


&nbsp;



# License
***
&nbsp;

This work is subject to the [MIT License](https://github.com/kevinzhongxu/CasOligo/LICENSE.txt).

&nbsp;
&nbsp;

&nbsp;

<hr />
<p style="text-align: center;">A work by <a href="https://kevinxuzhong.netlify.com/">Kevin Xu ZHONG</a></p>
<p style="text-align: center;"><span style="color: #808080;"><em>xzhong@eoas.ubc.ca</em></span></p>








