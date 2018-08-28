#!/usr/local/public/bin/Rscript --verbose
# version="1.1"

# date: 04-06-2013
# **Authors** Gildas Le Corguille  ABiMS - UPMC/CNRS - Station Biologique de Roscoff - gildas.lecorguille|at|sb-roscoff.fr 

# abims_hclust.r version 20130604

library(batch)
library(ctc)

hclust_metabolomics = function(file, method = "pearson", link = "ward", normalization=TRUE, keep.hclust=FALSE, sep=";", dec="."){

    if (sep=="tabulation") sep="\t"
    if (sep=="semicolon") sep=";"
    if (sep=="comma") sep=","

    # -- loading --
    data=read.table(file, header = TRUE, row.names=1, sep = sep, quote="\"", dec = dec,
		    fill = TRUE, comment.char="",na.strings = "NA")

    # -- Normalization: logratio --
    if (normalization) {
	    #meandata = apply(data,1,mean, na.rm=T)
	    #data = log2(data/meandata)
	    data=t(scale(t(data)))
        
        #AMAP: Unable to compute Hierarchical Clustering: missing values in distance matrix
        #Erreur dans hcluster(x, method = method, link = link) :
        #  Missing values in distance Matrix
        #Calls: do.call -> <Anonymous> -> hclust2treeview -> hcluster
        #Exécution arrêtée
        data[is.nan(data)] = 0
    }
    
    #Erreur dans `[.default`(xj, i) :
    #  les indices négatifs ne peuvent être mélangés qu'à des 0
    #Calls: do.call ... r2cdt -> [ -> [.data.frame -> [ -> [.factor -> NextMethod
    #Exécution arrêtée
    data = data[!apply(data,1,sum)==0,]

    # -- hclust / output files for TreeView --
    file="hclust.cdt"
    hclust2treeview(data,file=file, method = method, link = link, keep.hclust= keep.hclust)
    
    # -- output / return --
    system("zip -r hclust.zip hclust.*", ignore.stdout = TRUE)
}

listArguments = parseCommandArgs(evaluate=FALSE)
do.call(hclust_metabolomics, listArguments)
