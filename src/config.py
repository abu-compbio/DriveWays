
import time
from strategy.functions import *
# now = datetime.datetime.now()
# # String to save results on folder named with the date of the experiment
# today = '{}_{}_{}/'.format(now.day,now.month,now.year)


data_path = '../data/'
out_path = '../out/'
components = '../out/components/'
# condition = 'cov_over_mutex'
score = 'cov'
date = components+'04_11_2019'
# t = 1.4
# d = 4
N = 15
pathway_name = 'Kegg6.2'
# components_folder = '{}/{}_{}_t{}_{}_NoOverlapping/'.format(date,condition,score,t_threshold,pathway_name)
# hintedge_file = data_path+'hintedge_new.txt'
gene_patients_file =data_path+'gene_patient.txt'
k=3

mkdir(components)
mkdir(date)
# mkdir(components_folder)
mkdir(data_path)
