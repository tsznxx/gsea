#!/usr/bin/env python
# -- coding:utf-8 --
# Last-modified: 20 Jul 2017 05:19:34 PM
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
import pandas
import argparse

import gsea

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
    pr.add_argument("-c","--cls",dest="cls",type=str,metavar='sample.cls',required=True,help="Sample class file in .cls format.")
    pr.add_argument("-g","--gmt",dest="gmt",type=str,metavar="genesets.gmt",required=True,help="Gene set file in .gmt format.")
    pr.add_argument("-o",dest="outprefix",type=str,metavar="outprefix",required=True,help="GSEA output file prefix.")
    
    po = p.add_argument_group('Optional')
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

def GSEA(exprs,cls,gmt,prefix,ascending=False,method='ttest',weight=1,min_gene=15,max_gene=500,seed=1024,nperm=1000):
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
    grank = gsea.cal_rank(exprs,cls,prefix=prefix,ascending=ascending,method=method)
    gsets = gsea.read_genesets(gmt,min_gene=min_gene,max_gene=max_gene)
    gsea.calculateES3D(grank.index,gsets,prefix=prefix,ranking=grank['rnk'],weight=weight,nperm=nperm,seed=seed)

# def GSEAPreRank(genes,gsets,ascending=False,ranking=ranking,weight=0,min_gene=15,max_gene=500,seed=1024,nperm=1000)



# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main
# ------------------------------------

if __name__=="__main__":
    args = argParser()
    GSEA(args.exprs,args.cls,args.gmt,prefix=args.outprefix,ascending=args.ascend,method=args.method,weight=args.weight,min_gene=args.mingene,max_gene=args.maxgene,seed=args.seed,nperm=args.nperm)

