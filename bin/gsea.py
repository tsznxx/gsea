#!/usr/bin/env python
# -- coding:utf-8 --
# Last-modified: 21 Jul 2017 05:00:15 PM
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
import argparse

from scipy import stats
# from multiprocessing import Pool  

numpy.seterr(divide='ignore')

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

def argParser():
    ''' Parse arguments. '''
    p=argparse.ArgumentParser(description='Gene Set Enrichment Analysis. Implemented according to Subramanian, A., et al. (2005). "Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles." Proc Natl Acad Sci U S A 102(43): 15545-15550.',add_help=False,epilog='dependency numpy, scipy, pandas')

    pr = p.add_argument_group('Required')
    pr.add_argument("-e","--exprs",dest="exprs",type=str,metavar="gene_exprs.tsv", required=True, help="Gene expression matrix file in text format.")
    pr.add_argument("-g","--gmt",dest="gmt",type=str,metavar="genesets.gmt",required=True,help="Gene set file in .gmt format.")
    pr.add_argument("-o",dest="outprefix",type=str,metavar="outprefix",required=True,help="GSEA output file prefix.")

    po = p.add_argument_group('Optional')
    po.add_argument("-c","--cls",dest="cls",type=str,metavar='sample.cls',default=None,help="Sample class file in .clsformat.")
    po.add_argument("--lambda", dest="lamb",type=str,metavar="lamba x:x",default=None,help="A lambda expression to interpret the sample labels from the expression sample names.")
    po.add_argument("--method",dest="method",type=str,default='ttest',choices=['ttest'],help="Algorithms used for ranking metric calculation.")
    po.add_argument("--ascend",dest="ascend",action='store_true',default=False,help="Rank genes in ascending order. [Default is False.]")
    po.add_argument("--seed",dest="seed",type=int,metavar="1024",default=1024,help="Seed to generate random values.")
    po.add_argument("--min_gene",dest="mingene",type=int,metavar="15",default=15,help="Minimum number of genes in gene set.")
    po.add_argument("--max_gene",dest="maxgene",type=int,metavar="500",default=500,help="Maximum number of genes in gene set.")
    po.add_argument("-w",dest="weight",type=int,metavar="1",default=1,help="Weight of ranking metric.")
    po.add_argument("-n",dest='nperm',type=int,metavar='1000',default=1000,help="Number of permutations. [Default=1000]")
    if len(sys.argv)==1:
        sys.exit(p.print_help())
    args = p.parse_args()
    return args

def touchtime(lstr,ofh=sys.stderr):
    ofh.write("# {0}: {1}\n".format(time.ctime(),lstr))

def cal_rank(exprfile,clsfile,prefix,lamb=None,ascending=False,method='ttest'):
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
    if clsfile is None:
        cls = open(exprfile).next().rstrip().split('\t')[1:] # skip gene names
        if not lamb is None:
            fun = eval(lamb)
            cls = [fun(x) for x in cls]
        cls = numpy.array(cls)
    else:
        with open(clsfile) as fh:
            nsam, nlabel, _ = fh.next().split()
            _, ctl, trmt = fh.next().split()
            cls = numpy.array([int(i) for i in fh.next().split()])
    labels = list(set(cls))
    print labels
    ctl, trmt = cls[0], labels[1] if labels[0]==cls[0] else labels[0]
    touchtime("Labels identified: {0}:{1}, {2}:{3} ...".format(ctl,(cls==ctl).sum(),trmt,(cls==trmt).sum()))
    
    # ranking
    touchtime("Calculate ranking metric using {0} method ...".format(method))
    idx = cls == cls[0] # first encountered label as control
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
        touchtime("Weighted correclation vector is suppressed because either ranking metric is not provided or weight is set to 0.")
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
    rdf = rdf.sort_values('NES',ascending=False)
    rdf.to_csv("{0}_GSEA.tsv".format(prefix),sep="\t")
    return rdf


#def _cal_es_pval(ind,rs=numpy.random.RandomState(),niter=1000):
#    '''
#    Private function for permutation and pvalue calculation using empirical CDF.
#    Parameters:
#        ind: pandas.Series of bools
#            indicator of gene in geneset
#        rs: numpy.random.RandomState()
#            randome state generator
#        niter: int
#            number of iteration
#    Returns:
#        gsid: string
#            gene set id. Used to keep track of the gene sets
#        es: float
#            enrichment score
#        pval: float
#            pvalue inferred from empirical CDF.
#    NOTE: 
#        we do some tricks here.
#        1. Assign hits and nonhits with ratios by: indmat*(nMiss_ratio+nHit_ratio)-nMiss_ratio
#        2. indmat is N by niter+1. The last element of indmat is not shuffled.
#    '''
#
#    # assign ratio to hits and nonhits
#    N = len(ind)
#    nHit = sum(ind) # Number of hits
#    if nHit == 0:
#        return 0.,1.
#    nHit_ratio= 1./nHit # ratio of hits, sum(hits)=1
#    nMiss_ratio = 1./(N-nHit) # ratio of misses, sum(nonhits)=-1
#
#    # generate permutation matrix
#    indmat = ind.repeat(niter+1).reshape(N,niter+1).T
#    numpy.apply_along_axis(rs.shuffle ,1,indmat[:-1]) # the last one is not shuffled.
#
#    # calculate ens
#    es = (ind*(nMiss_ratio+nHit_ratio)-nMiss_ratio).cumsum() # 
#    indmat = indmat*(nMiss_ratio+nHit_ratio)-nMiss_ratio
#    ens = numpy.cumsum(indmat,1)
#    EMAX, EMIN = numpy.max(ens,1), numpy.min(ens,1)
#    esl = numpy.where(numpy.abs(EMAX)>numpy.abs(EMIN),EMAX,EMIN)
#    es, esl = esl[-1], esl[:-1]
#    
#    # Calculate pvalue using empirical CDF (last element is the enrichment score to test) 
#    # Normalize the ES score by average ES
#    # possible errors: divided by zero or nan. Set pval=1.0 and nes=0.0
#    try:
#        if es >=0:
#            pval = float(numpy.sum(esl>es))/numpy.sum(esl>=0)
#            nes = es/esl[esl>=0].mean()
#        else:
#            pval = float(numpy.sum(esl<es))/numpy.sum(esl<0)
#            nes = -es/esl[esl<0].mean() 
#    except:
#        pval, nes = 1.0, 0.0
#    return es,nes,pval
#def calculateES(grank,gsetfile,nproc=10,niter=1000,seed=1024,min_gene=15):
#    '''
#    Do permuation on genesets and calculate pvalue using empirical CDF.
#    Parameters:
#        grank: pandas.DataFrame
#            gene expression data frame sorted by 'ranking' socres.
#        gsetfile: string
#            gene set file. 1st column is the gene set ID. 2nd column is the gene set description and the rest is the genes in the gene set.
#
#    '''
#    # get geneset hit matrix
#    gsets = {}
#    descriptions = {}
#    for line in open(gsetfile):
#        items = line.rstrip().split('\t')
#        descriptions[items[0]] = items[1]
#        if len(items)>= min_gene:
#            gsets[items[0]] = items[2:]
#    rdf = pandas.DataFrame({'description':descriptions.values()},index=descriptions.keys())
#    rdf = rdf.sort_index()
#    rdf.name = "gset"
#
#    # run in multiple processes
#    rs = numpy.random.RandomState(seed)
#    pool = Pool(processes=nproc) 
#    ts = []
#    for gsid in rdf.index[:]:
#        ts.append(pool.apply_async(_cal_es_pval,args=(grank.isin(gsets[gsid]),rs,niter)))
#    pool.close()
#    pool.join()
#
#    # summarize results to dataframe.
#    rdf['ES'], rdf['NES'],rdf['pvalue'] = zip(*[t.get() for t in ts])
#    rdf = rdf.sort_values('NES',ascending=False)
#    return rdf

def GSEA(exprs,cls,gmt,prefix,ascending=False,method='ttest',lamb=None,weight=1,min_gene=15,max_gene=500,seed=1024,nperm=1000):
    '''
    Run GSEA algorithm.
    Parameters:
        exprs: pandas.DataFrame
            gene expression matrrix in gene by sample format.
        cls: list or pandas.Series
            sample class labels.
        gmt: string
            gene set matrix file.
        ascending: bool
            Rank genes in ascending order. [Default is False.]
        method: string
            method used for ranking calculation.
        lamb: string
            lambda expression to interpret class labels from sample names.
        weight: float
            GSEA's weighted score. weighted_r = r**weight. This is to balance the weight of the ranking metric.
        min_gene, max_gene: int
            minimum/maximum number of genes required in the get set.
        seed: int
            used in numpy.randome.RandomState(seed).
        nperm: int
            number of permutations.
    Returns:
        rdf: pandas.DataFrame
            result dataframe contains: gsetid (as index), description, ES, NES, pvalue, FDR
    '''
    grank = cal_rank(exprs,cls,prefix=prefix,ascending=ascending,lamb=lamb,method=method)
    gsets = read_genesets(gmt,min_gene=min_gene,max_gene=max_gene)
    calculateES3D(grank.index,gsets,prefix=prefix,ranking=grank['rnk'],weight=weight,nperm=nperm,seed=seed)


# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main
# ------------------------------------

if __name__=="__main__":
    args = argParser()
    GSEA(args.exprs,args.cls,args.gmt,prefix=args.outprefix,ascending=args.ascend,lamb=args.lamb,method=args.method,weight=args.weight,min_gene=args.mingene,max_gene=args.maxgene,seed=args.seed,nperm=args.nperm)

