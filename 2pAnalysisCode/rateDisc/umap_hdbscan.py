# -*- coding: utf-8 -*-
"""
Created on Sun May 10 15:42:38 2020

@author: Simon
"""

import numpy as np
import h5py
import umap
import hdbscan
import scipy.io as sio


allPath = ['Q:\\BpodImager\\umapClust\\20200510T1616_umap_area2.mat',
           'Q:\\BpodImager\\umapClust\\20200510T1616_umap_area3.mat',
           'Q:\\BpodImager\\umapClust\\20200510T1616_umap_area4.mat',
           'Q:\\BpodImager\\umapClust\\20200510T1616_umap_area5.mat',
           'Q:\\BpodImager\\umapClust\\20200510T1616_umap_area6.mat',
           'Q:\\BpodImager\\umapClust\\20200510T1616_umap_area7.mat',
           'Q:\\BpodImager\\umapClust\\20200510T1616_umap_area8.mat',
           'Q:\\BpodImager\\umapClust\\20200510T1616_umap_area9.mat',
           'Q:\\BpodImager\\umapClust\\20200510T1616_umap_area10.mat',
           'Q:\\BpodImager\\umapClust\\20200510T1616_umap_area246.mat',
           'Q:\\BpodImager\\umapClust\\20200510T1616_umap_area247.mat',
           'Q:\\BpodImager\\umapClust\\20200510T1616_umap_area248.mat',
           'Q:\\BpodImager\\umapClust\\20200510T1616_umap_area249.mat',
           'Q:\\BpodImager\\umapClust\\20200510T1616_umap_area250.mat',
           'Q:\\BpodImager\\umapClust\\20200510T1616_umap_area251.mat',
           'Q:\\BpodImager\\umapClust\\20200510T1616_umap_area252.mat',
           'Q:\\BpodImager\\umapClust\\20200510T1616_umap_area253.mat',
           'Q:\\BpodImager\\umapClust\\20200510T1616_umap_area254.mat']

# allPath = ['Q:\\BpodImager\\umapClust\\20200512T0009_umap_all.mat']


for cPath in allPath:
    # load data
    print(cPath)
    g = h5py.File(cPath,'r')
    X = np.array(g['X'])
    g.close
    
    # run umap twice, once for visualization and one in higherD for clustering
    reducer = umap.UMAP()
    umapOut = reducer.fit_transform(X)
    print('Umap vals: ' + str(umapOut.shape))
    
    umapClusterable = umap.UMAP(
        n_neighbors=30,
        min_dist=0.0,
        n_components=15,
        random_state=42,
        ).fit_transform(X)
    
    # run hdbscan on highD embedding
    hdbLabels = hdbscan.HDBSCAN(min_samples=15,min_cluster_size=15).fit_predict(umapClusterable)
    
    sio.savemat(cPath.replace('.mat','_clustered.mat'),
                    {'umapOut':umapOut,
                     'umapClusterable':umapClusterable,
                     'hdbLabels':hdbLabels,
                     })