#!/usr/bin/env python
# -- coding:utf-8 --
# Last-modified: 20 Jul 2017 05:25:57 PM
#
#         Module/Scripts Description
# 
# Copyright (c) 2016 Yunfei Wang <yfwang0405@gmail.com>
# 
# This code is free software; you can redistribute it and/or modify it
# under the terms of the BSD License (see the file COPYING included with
# the distribution).
# 
# @status:  experimental
# @version: 1.0.0
# @author:  Yunfei Wang
# @contact: yfwang0405@gmail.com

# ------------------------------------
# python modules
# ------------------------------------

import os
import sys
import time
import numpy
import pandas
from scipy import stats
# from multiprocessing import Pool  

numpy.seterr(divide='ignore')

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

def touchtime(lstr,ofh=sys.stderr):
    ofh.write("# {0}: {1}\n".format(time.ctime(),lstr))

def cal_rank(exprfile,clsfile,prefix,ascending=False,method='ttest'):
    '''
    Calculate rank from gene expression file, class file.
    Parameters:
        exprs: pandas.DataFrame
            gene expression matrrix in gene by sample format.
        cls: list or pandas.Series
            sample class labels. Categorical (e.g tumor vs normal) class file format (*.cls).
            http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_forma
            ts#CLS:_Categorical_.28e.g_tumor_vs_normal.29_class_file_format_.28.2A.cls.29
        ascending: bool
            False: up-regulated genes rank first.
        method: string
          Algorithms used for ranking metric calculation. Choices from ['ttest',...]
    '''
    # sample class
    touchtime("Read gene expression file ...")
    # gene expression
    expr = pandas.read_table(exprfile,index_col=0)
    
    touchtime("Read sample class file ...")
    with open(clsfile) as fh:
        nsam, nlabel, _ = fh.next().split()
        _, ctl, trmt = fh.next().split()
        cls = numpy.array([int(i) for i in fh.next().split()])
    
    # ranking
    touchtime("Calculate ranking metric using {0} method ...".format(method))
    idx = cls == pandas.unique(cls)[0]
    grank = pandas.DataFrame(index=expr.index)
    if method == 'ttest':
        tscores, pvalues = stats.ttest_ind(expr.loc[:,idx],expr.loc[:,~idx],axis=1)
        grank['tscores'] = tscores
        grank['rnk'] = -numpy.log10(pvalues)*numpy.sign(tscores)
    else:
        pass # reserved for future use
    grank = grank.sort_values('rnk',ascending=ascending)
        
    touchtime("Save ranking metric to file: {0}_ranking_metric.tsv ...".format(prefix))
    grank.to_csv("{0}_ranking_metric.tsv".format(prefix),sep='\t')
    return grank

def read_genesets(gmtfile,min_gene=15,max_gene=500):
    '''
    Read Geneset Matrix Table (.gmt) into dict. 
    Format definition:
      http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29)
    Parameters:
        gmtfile: string
            geneset in gmt format. 1st column is the gene set ID. 2nd column is the gene set description and the rest are the genes in the gene set.
        min_gene, max_gene: int
            minimum/maximum number of genes required in the get set.            
    Returns:
        gsets: dict
            gene set dict. gsets = {'gsid':['description',[genes]]}
    '''
    # read .gmt file
    touchtime("Reading geneset matrix file ...")
    gsets = {}
    for line in open(gmtfile):
        items = line.rstrip().split('\t')
        if len(items)>= min_gene:
            gsets[items[0]] = [items[1],items[2:]]
    return gsets

def calculateES3D(grank,gsets,prefix,ranking=None,weight=1,nperm=1000,seed=1024):
    '''
    Do permuation on genesets and calculate pvalue using empirical CDF.
    Parameters:
        grank: pandas.DataFrame
            index is the sorted genes. The 'ranking' column contains the ranking metric.
        gsets: dict
            gene set dict. gsets = {'gsid':['description',[genes]]}
        prefix: string
            prefix of result file.    
        ranking: None or a list
            None: no correlation vector.
            list: the correlation vector is provided as a list.
        weight: float
            GSEA's weighted score. weighted_r = r**weight. This is to balance the weight of the ranking metric.
        seed: int
            used in numpy.randome.RandomState(seed).
        nperm: int
            number of permutations.
    Returns:
        rdf: pandas.DataFrame
            GSEA results with geneset information.
    '''
    # Generate result DataFrame
    touchtime("Generate result data frame ...")
    rdf = pandas.DataFrame({'description':[val[0] for val in gsets.values()]},index=gsets.keys())
    rdf = rdf.sort_index()
    rdf.name = "gset"
    
    # weighted correlation vector (normalized to avoid non-positive values)
    touchtime("Calculate weighted correlation vector ...")
    # correlation vector
    N, M = grank.shape[0], rdf.shape[0] 
    if ranking is None or weight==0:
        cor_vec = numpy.repeat(1.0,N)
    else: # list
        cor_vec = numpy.abs(ranking**weight)  # correlation vector

    # generate inhit matrix, permutation matrix
    touchtime("Create permutation matrix ...")
    hit_mat = numpy.array([grank.isin(gsets[gsid][1]) for gsid in rdf.index]) # M genestes by N genes
    perm_mat = numpy.repeat(hit_mat,nperm+1).reshape((M,N,nperm+1))    
    
    touchtime("Shuffling permutation matrix ...")
    # shuffle, last matrix is not shuffled
    rs = numpy.random.RandomState(seed)    
    numpy.apply_along_axis(lambda x: numpy.apply_along_axis(rs.shuffle,0,x),1,perm_mat[:,:,:-1]) 
        
    # assign ratio to hits and nonhits   
    touchtime("Calculate step factors ...")
    miss_perm_mat = perm_mat==0
    perm_mat = perm_mat* cor_vec[numpy.newaxis,:,numpy.newaxis]    # multiply by correlation vector
    hit_factor = 1./perm_mat.sum(1) # ratio of hits, sum(hits)=1 
    hit_factor[numpy.isinf(hit_factor)] = 0 # set inf values to zero. if cor_vec is scaled to (0,1], there should be no inf values.
    miss_factor = 1./miss_perm_mat.sum(1) # ratio of misses, sum(nonhits)=-1
    
    touchtime("Calculate enrichment scores ...")
    scores = (perm_mat*hit_factor[:,numpy.newaxis,:] - miss_perm_mat*miss_factor[:,numpy.newaxis,:]).cumsum(1)
    emax,emin = numpy.max(scores,1), numpy.min(scores,1)
    esall = numpy.where(numpy.abs(emax)>numpy.abs(emin),emax,emin)
    es, esnull = esall[:,-1], esall[:,:-1]
    
    touchtime("Calculate p values using empirical CDF ...")
    sumPos, sumNeg = (esnull>=0).sum(1).astype(float), (esnull<0).sum(1).astype(float)
    pvalPos = (esnull>=es[:,numpy.newaxis]).sum(1)/sumPos
    pvalNeg = (esnull< es[:,numpy.newaxis]).sum(1)/sumNeg
    pval = numpy.where(es>=0,pvalPos,pvalNeg)    
    
    touchtime("Normalize enrichment scores ...")
    avgPos = (esnull*(esnull>=0)).sum(1)/sumPos
    avgNeg = (esnull*(esnull<0)).sum(1)/sumNeg
    norm_es  = numpy.where(es>=0,es/avgPos,-es/avgNeg)
    norm_esnull = numpy.where(esnull>=0,esnull/avgPos[:,numpy.newaxis],-esnull/avgNeg[:,numpy.newaxis])
    
    touchtime("Calculate FDR ...")
    # percentage calculattion
    #     norm_es >= 0: (N-idx_in_list)/nPos
    #     norm_es <  0: idx_in_list/nNeg
    # observed percentage
    sorted_nes = numpy.sort(norm_es)
    nNeg = numpy.searchsorted(sorted_nes,0).astype(float)
    nPos = sorted_nes.shape[0] - nNeg
    norm_es_idx = numpy.searchsorted(sorted_nes,norm_es)
    obs_percent = numpy.where(norm_es>=0,(sorted_nes.shape[0]-norm_es_idx)/nPos,norm_es_idx/nNeg)    
    
    # expected percentage
    norm_esnull = numpy.sort(norm_esnull.flatten())
    nNeg = numpy.searchsorted(norm_esnull,0).astype(float)
    nPos = norm_esnull.shape[0] - nNeg
    norm_es_idx = numpy.searchsorted(norm_esnull,norm_es)
    expect_percent = numpy.where(norm_es>=0,(norm_esnull.shape[0]-norm_es_idx)/nPos,norm_es_idx/nNeg)
    
    # FDR
    fdr = expect_percent/obs_percent
    fdr[fdr>1.] = 1.0    
    fdr[expect_percent==0] = 0. # 
    
    # summarize results to dataframe.
    touchtime("Summarize results ...")
    rdf['ES'], rdf['NES'],rdf['pvalue'], rdf['FDR'] = es, norm_es, pval, fdr
    ngenes, gstrs = zip(*[(len(gsets[gsid][1]),",".join(gsets[gsid][1])) for gsid in rdf.index])
    rdf['nGenes'], rdf['nHits'], rdf['genes'] = ngenes, hit_mat.sum(1), gstrs
    
    # save results
    touchtime("Save GSEA results to {0}_GSEA.tsv ...".format(prefix))
    rdf.to_csv("{0}_GSEA.tsv".format(prefix),sep="\t")
    return rdf

# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main
# ------------------------------------

if __name__=="__main__":
    if len(sys.argv)==1:
        sys.exit("Example:"+sys.argv[0]+" infile1 infile2 [infile3 ...]\n{0} ".format(DESCRIPTION))

