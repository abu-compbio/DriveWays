from utils import *
import glob
import os
import math
import operator
import seaborn as sns
sns.set(style="whitegrid")
from tqdm import trange, tqdm
# from matplotlib import pyplot as plt
#from matplotlib import rcParams
# from pylab import *
# import matplotlib as mpl
# mpl.style.use('seaborn')
#plt.figure(figsize=(25, 20))
def Logpval_per_comp():
    models = [
            'cov_seed_C_cov_over_mutex_S_cov','cov_seed_C_cov_over_mutex_S_cov*mutex',
            'cov_seed_C_cov_2_S_cov','cov_seed_C_cov_2_S_cov*mutex',
            'cov_seed_C_cov*mutex_S_cov', 'cov_seed_C_cov*mutex_S_cov*mutex',

            'cov*mutex_seed_C_cov_over_mutex_S_cov','cov*mutex_seed_C_cov_over_mutex_S_cov*mutex',
            'cov*mutex_seed_C_cov_2_S_cov', 'cov*mutex_seed_C_cov_2_S_cov*mutex',
            'cov*mutex_seed_C_cov*mutex_S_cov','cov*mutex_seed_C_cov*mutex_S_cov*mutex',
            ]

    # models= ['mexcowalk','cov_over_mutex_cov_0.1']


    labels = models
    model2label = {m:l for m,l in zip(models,labels )}

    model2w = {}
    with open('../out/components/11_03_2019/cancer_subtype_test_results/Logpval_per_comp.txt') as f:
        lines = f.readlines()
        Ns = [int(n) for n in lines[0].rstrip().split('\t')]
        for line in lines[1:]:
            line = line.rstrip().split('\t')
            models_ = line[0]
            W = [float(w) for w in line[1:] if w != 'N/A']
            model2w[models_] = W


    axes(frameon=0)
    legend_ = []
    colors = [None]*10+['r','b']
    for m,c in zip(models,colors):
        legend_.append(model2label[m])
        try:
            plot(Ns, model2w[m], '-o', color = c)
        except:
            print(m, len(Ns), len(model2w[m]))

    art = []

    axis_font = {'fontname':'Times New Roman', 'size':'15'}

    legend = plt.legend(legend_, loc=8,fancybox=True, fontsize= 'medium', framealpha=0,
                        edgecolor = 'b', ncol= 4, bbox_to_anchor=(0.5,-0.25))
    art.append(legend)
    frame = legend.get_frame()

    plt.xlabel('total_genes', **axis_font)
    plt.ylabel('Cancer Type Specificity Score (CTSS)', **axis_font)

    xtick = list(range(100,2600,200))#+[554]+ list(range(600,800, 200))+[806]+list(range(900, 2600, 200))
    xticks(xtick, fontsize='x-small', **axis_font)
    # plt.ylim(-0.01, 30)
    # axes(frameon=0)
    grid()
    plt.savefig('../out/components/11_03_2019/cancer_subtype_test_results/Logpval.pdf',format = 'pdf', additional_artists=art,
                bbox_inches="tight", dpi = 800)
    plt.close()


def max_size():


    models  = ["hotnet2",
    "memcover_v1","memcover_v2","memcover_v3",
    "mutex_t07_nsep_cov",
    'hier_hotnet'
    ]

    legend_ = [
        "Hotnet2",
        "MEMCover_v1",
        "MEMCover_v2",
        "MEMCover_v3",
        "MEXCOwalk",
        "HierHotnet_v1",
        "HierHotnet_v2",
        ]
    component_path = '../hint/out/connected_components_isolarge_n2500_whh/'
    N = list(range(100,600,100))+[554]+ list(range(600,800, 100))+[806]+list(range(900, 2600, 100))
    memcover_v3_range = list(range(100,600,100))+[554]+ list(range(600,800, 100))+[806]+list(range(900, 1700, 100))
    for model in models:
        delimiter = '\t' if 'hier' in model or  'memcover' in model else ' '
        res = []
        max_size = []
        for n in N:
            if (n not in [554, 806] and model == 'hier_hotnet') or (n not in memcover_v3_range and model == 'memcover_v3'):
                continue
            else:
                file = glob.glob(component_path +'{}/cc_n{}_*'.format(model, n))[0]
                # print('file')
                with open(file, 'r') as f:
                    lines = [len(l.strip().split(delimiter)) for l in f.readlines()]
                    max_ = max(lines)
                    max_modules = [1 for m in lines if m == max_]

                    res.append(max_*len(max_modules)/n)
                    max_size.append(max_)
        print(model,' ', max_size)
        axes(frameon=0)
        if model == 'hier_hotnet':
            plot([806], res[1], 'k*', markersize=12)
            plot([554], res[0], 'C8*', markersize=12)
        elif model == 'memcover_v3':
            plot(memcover_v3_range, res, '-o')
        else:
            plot(N, res,  '-o')

    art = []


    legend = plt.legend(legend_, loc=8,fancybox=True, fontsize= 'small', framealpha=0,
                        edgecolor = 'b', ncol= 4, bbox_to_anchor=(0.5,-0.3))
    art.append(legend)
    frame = legend.get_frame()
    plt.xlabel('total_genes')
    plt.ylabel('Percentage of Genes in Max Sized Modules')

    xtick = list(range(100,2600,200))#+[554]+ list(range(600,800, 200))+[806]+list(range(900, 2600, 200))
    xticks(xtick, fontsize='x-small')
    # plt.ylim(-0.01, 40)
    # axes(frameon=0)
    grid()
    plt.savefig('../results/main plots/max_size_over_N.pdf',format = 'pdf', additional_artists=art,
                bbox_inches="tight", dpi = 800)
    plt.close()



def min_size_modules():


    models  = ["hotnet2",
    "memcover_v1","memcover_v2","memcover_v3",
    "mutex_t07_nsep_cov",
    'hier_hotnet'
    ]

    legend_ = [
        "Hotnet2",
        "MEMCover_v1",
        "MEMCover_v2",
        "MEMCover_v3",
        "MEXCOwalk",
        "HierHotnet_v1",
        "HierHotnet_v2",
        ]
    # print(glob.glob('../results/results/*'))
    component_path = '../hint/out/connected_components_isolarge_n2500_whh/'
    N = list(range(100,600,100))+[554]+ list(range(600,800, 100))+[806]+list(range(900, 2600, 100))
    memcover_v3_range = list(range(100,600,100))+[554]+ list(range(600,800, 100))+[806]+list(range(900, 1700, 100))
    for model in models:
        delimiter = '\t' if 'hier' in model or  'memcover' in model else ' '
        res = []
        min_size = []
        for n in N:
            if (n not in [554, 806] and model == 'hier_hotnet') or (n not in memcover_v3_range and model == 'memcover_v3'):
                continue
            else:
                file = glob.glob(component_path +'{}/cc_n{}_*'.format(model, n))[0]
                # print('file')
                with open(file, 'r') as f:
                    lines = [len(l.strip().split(delimiter)) for l in f.readlines()]
                    min_ = min(lines)
                    min_modules = [1 for m in lines if m == min_]

                    res.append(min_*len(min_modules)/n)
                    min_size.append(min_)

        print(model,' ', min_size)
        # print(len(min_size), ' ', len(N))
        axes(frameon=0)
        if model == 'hier_hotnet':
            plot([806], res[1], 'k*', markersize=12)
            plot([554], res[0], 'C8*', markersize=12)
        elif model == 'memcover_v3':
            plot(memcover_v3_range, res, '-o')
        else:
            plot(N, res,  '-o')

    art = []


    legend = plt.legend(legend_, loc=8,fancybox=True, fontsize= 'small', framealpha=0,
                        edgecolor = 'b', ncol= 4, bbox_to_anchor=(0.5,-0.3))
    art.append(legend)
    frame = legend.get_frame()
    plt.xlabel('total_genes')
    plt.ylabel('Percentage of Genes in Min Sized Modules')

    xtick = list(range(100,2600,200))#+[554]+ list(range(600,800, 200))+[806]+list(range(900, 2600, 200))
    xticks(xtick, fontsize='x-small')
    # plt.ylim(-0.01, 40)
    # axes(frameon=0)
    grid()
    plt.savefig('../results/main plots/min_size_over_N.pdf',format = 'pdf', additional_artists=art,
                bbox_inches="tight", dpi = 800)
    plt.close()
def mmr_plot():
    names = ['ClusterOne', 'MexCoGrowth', 'Hotnet2', 'MEMCover', 'MexCoWalk']
    kegg_mmr = [0.02862,0.0850,0.0312,0.05271,0.00909]
    Reactome_mmr = [0.0428,0.0923,0.0068,0.06251,0.0388]
    biocarta_mmr = [0.03340,0.1331,0.0494,0.0691,0.0376]

    kegg_ratio = np.array([[105,102,136,155,171],[1291,602,1770,1386,1770]])
    reactome_ratio = np.array([[145,145,214,216,246],[2212,1199,3368,2649,3368]])
    biocata_ratio = np.array([[62,60,67,90,91],[893,439,1173,921,1173]])

    kegg_ratio = kegg_ratio[0]/kegg_ratio[1]
    reactome_ratio = reactome_ratio[0]/reactome_ratio[1]
    biocata_ratio = biocata_ratio[0]/biocata_ratio[1]

    for a in [kegg_mmr,Reactome_mmr,biocarta_mmr]:
        sns.barplot(names, a)


mmr_plot()

# diff_mutex_thresh_plot()

# Logpval_per_comp()

# print('Max size for each N')
# max_size()
# print('Min size for each N')
# min_size_modules()
