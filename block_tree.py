import queue
import os
import re
from collections import deque


def init_newick(file_txt): 
    
    list1  = [[r':0.\d+|:1', ''], [r'\s*,\s*', ','],  [r'\s*\(\s*', '('], [r'\s*\)\s*', ')']]
    result = file_txt
    for pat in list1:
        result = re.sub(pat[0],pat[1],result)
    result = result.strip()
    if result[-1]!=')':
        result = result[:-1]
   
    return result
  
    

def trans_newick(tree_dict, root, tree_taxa): 
       
    if root[0] != root[1]:
        if tree_taxa[root[0]][1] > 0:
            tree_taxa[root[0]][1] -= 1
        if tree_taxa[root[1]][2] > 0:
            tree_taxa[root[1]][2] -= 1            
        child = []
        left = 0
        right = 0
        start = root[0] 
        
        for i in range(root[0], root[1] + 1):            
            left += tree_taxa[i][1]
            right += tree_taxa[i][2]
            if left == right:
                child.append([start,i]) 
                start = i + 1        
                  
        len1 = len(child)          
        for i in range(len1):
            lr = (child[i][0], child[i][1])         
            tree_dict[root][-1].append(lr)
            a = [tree_taxa[j][0] for j in range(child[i][0], child[i][1]+1)]
            a = frozenset(a)            
            tree_dict[lr] = [a, root, []]              
            trans_newick(tree_dict, lr, tree_taxa)
      
         
        
def read_files(cur_dir):     
    
    file_in = os.listdir(cur_dir)
    file_num = len(file_in)         
    tree_taxa_dict = {i:[] for i in range(file_num)}
    tree_dict_list = [dict() for i in range(file_num)]
    tree_file_list = [[] for i in range(file_num)]
    origin_to_taxa = dict() 
    taxa_to_origin = dict()
    
    for i in range(file_num):
        file = cur_dir + file_in[i]
        tree_file_list[i] = file_in[i]
        with open(file, 'r') as f:
            file_txt = f.read()
        a = init_newick(file_txt)   
        a = a.split(',')
        len1 = len(a)
        tree_taxa_dict[i] = [[] for j in range(len1)]
        for j  in range(len1):
            taxa = re.sub(r'\(|\)', '', a[j])
            origin_to_taxa[taxa] = -1
            tree_taxa_dict[i][j] = [taxa, a[j].count('('), a[j].count(')')] 
      
    idx = 1       
    for i in origin_to_taxa.keys():
        origin_to_taxa[i] = idx 
        taxa_to_origin[idx] = i
        idx += 1
        
    for i in range(file_num):       
        for j in range(len(tree_taxa_dict[i])): 
            tree_taxa_dict[i][j][0] = origin_to_taxa[tree_taxa_dict[i][j][0]] 
        taxa_set = [k[0] for k in tree_taxa_dict[i]] 
        root = (0, len(tree_taxa_dict[i])-1)        
        tree_dict_list[i][root] = [frozenset(taxa_set), -1, []]            
        trans_newick(tree_dict_list[i], root, tree_taxa_dict[i])
        tree_taxa_dict[i] = taxa_set
        
    return   tree_file_list, tree_dict_list, tree_taxa_dict, taxa_to_origin
        
        

def init_file_cluster(tree_dict_list, tree_taxa_dict, tree_test, taxa_occur): 
 
    cluster_in = [ ]
    common_taxa = [ ]       
    a = []
    for i in tree_test:
        a.extend(tree_taxa_dict[i])
    a1 = list(set(a))
    for i in a1:
        if a.count(i) >= taxa_occur:
            common_taxa.append(i)            

    for i in tree_test:                
        for j in tree_dict_list[i].keys():
            k = tree_dict_list[i][j][0].intersection(frozenset(common_taxa))
            if len(k) > 1 and len(k) < len(common_taxa):
                cluster_in.append(k)
            
    return cluster_in, len(common_taxa), len(cluster_in)



def connected_graph(conflict_dict):           
 
    num = len(conflict_dict)
    visit = [0 for u in range(num)]   
    component = [[]]  
    idx = 0 
    for i in range(num):
        if visit[i] == 0:
            if len(conflict_dict[i]) == 0:
                component[0].append(i)
            else: 
                idx += 1
                q = queue.Queue()
                q.put(i)
                visit[i] = 1  
                component.append([i])    
                while not q.empty():
                    a = q.get()  
                    for j in conflict_dict[a]:
                        if visit[j] == 0:
                            q.put(j)
                            visit[j] = 1
                            component[idx].append(j)    
               
    return  component           
              

    
def get_component(cluster_in):      
  
    num = len(cluster_in)    
    conflict_dict = {i:[] for i in range(num)}   
    
    for i in range(num - 1):
        for j in range(i + 1, num):
            len1 = len(cluster_in[i].intersection(cluster_in[j]))
            if len1 > 0 and len1 < len(cluster_in[i]) and len1 < len(cluster_in[j]):
                conflict_dict[i].append(j)
                conflict_dict[j].append(i)                     
                
    component = connected_graph(conflict_dict)      
    return component



def get_component_stat(cluster_comp):     
    
    num = len(cluster_comp)    
    stat = []
    for i in range(num - 1):
        for j in range(i + 1, num):
            s = cluster_comp[i].intersection(cluster_comp[j])
            s1 = cluster_comp[i].difference(s)
            s2 = cluster_comp[j].difference(s)
            if len(s) > 0 and len(s1) > 0 and len(s2) > 0:
                stat.append(frozenset([s, s1, s2]))
    
    stat = list(frozenset(stat))   
    for i in range(len(stat)):
        stat[i] = list(stat[i])                     
  
    return stat


  
def super_taxa_component(cluster_list):
    
    taxa_super = dict()
    cluster_comp = []   
    taxa_comp = []    
    len1 = len(cluster_list)
    
    a = []
    for i in cluster_list:
        a.extend(i)
    a = list(set(a))     
    taxon_occur = {i:[] for i in a}      
    for idx in range(len1):  
        for i in cluster_list[idx]:
            taxon_occur[i].append(idx)    
     
    temp = dict()     
    for key, value in taxon_occur.items(): 
        a = frozenset(value)
        if not (a in temp):              
            temp[a] = [key]     
        else: 
            temp[a].append(key)    
    for value in temp.values():  
        a = min(value)
        taxa_super[a] = value
        taxa_comp.append(a)
   
    for idx in range(len1):    
        a = cluster_list[idx].intersection(frozenset(taxa_comp))
        if len(a) > 1 and len(a) < len(taxa_comp):
            cluster_comp.append(a)    
                     
    return taxa_super, taxa_comp, cluster_comp




def component_summary(cluster_in):  
    
    component = get_component(cluster_in) 
    num = len(component)
    taxa_super_list = [[] for i in range(num)]
    taxa_comp_list = [[] for i in range(num)]
    cluster_comp_list = [[] for i in range(num)]      
    state_dict = [[] for i in range(num)]
 
    for idx in range(1,num):
        cluster_comp = [cluster_in[i] for i in component[idx]]             
        taxa_super_list[idx], taxa_comp_list[idx], cluster_comp_list[idx] = super_taxa_component(cluster_comp)
        state_dict[idx] = get_component_stat(cluster_comp_list[idx])   
  
    for i in component[0]:
        cluster_comp_list[0].append(cluster_in[i])
        taxa_comp_list[0].extend(cluster_in[i])
    taxa_comp_list[0] = list(set(taxa_comp_list[0])) 
    
    return  taxa_super_list, taxa_comp_list, cluster_comp_list, state_dict



def block_dict_insert(idx_block, block_idx, block, block_maxid):
   
    idx = block_maxid + 1
    if block in block_idx:
        idx = block_idx[block]
        idx_block[idx][1] += 1
    else:
        block_idx[block] = idx
        idx_block[idx] = [block,1]
        block_maxid = idx
    
    return idx, block_maxid



def block_dict_del(idx_block, block_idx, idx):
    
    if idx_block[idx][1] <= 1:
        a = idx_block[idx][0]
        del block_idx[a]
        del idx_block[idx]            
    else:
        idx_block[idx][1] -= 1
            



def block_union_check(taxa_set, block_set, set_len, idx_block):
    
    overlap = [1 for i in range(set_len)] 
    sum1 = []
    tf = False 
    lap_tf = 0 
    
    for i in range(set_len):  
        a = idx_block[block_set[i]][0]
        k = taxa_set.intersection(a)         
        if len(k) > 0:            
            sum1.extend(k)
            if len(k) < len(a):
                overlap[i] = k
                lap_tf = 1
    
    if lap_tf == 0 and len(sum1) == len(taxa_set):
        tf = True
                   
    return tf, overlap, sum1

    

def block_resolve(incomp_stat, block_set, idx_block, block_idx, block_maxid):
    
    overlap = [0, 0, 0] 
    block_sum = [0, 0, 0]
    set_len = len(block_set)
    flag = False
    for i in range(0, 3):         
        flag, overlap[i], block_sum[i] = block_union_check(incomp_stat[i], block_set, set_len, idx_block)
        if flag == True:
            break     
      
    blocks = 0    
    if flag == False:  
        blocks = [[], [], []] 
        for i in range(0, 3): 
            blocks[i].extend(block_set) 
            for bi in range(set_len):
                if overlap[i][bi] == 1:
                    idx_block[block_set[bi]][1] += 1
                else:
                    a = blocks[i][bi]
                    idx, block_maxid = block_dict_insert(idx_block, block_idx, overlap[i][bi], block_maxid)
                    blocks[i][bi] = idx
                
                    a = idx_block[block_set[bi]][0].difference(overlap[i][bi])
                    idx, block_maxid = block_dict_insert(idx_block, block_idx, a, block_maxid)
                    blocks[i].append(idx)        
                    
            if len(block_sum[i]) < len(incomp_stat[i]):   
                a = incomp_stat[i].difference(frozenset(block_sum[i]))   
                idx, block_maxid = block_dict_insert(idx_block, block_idx, a, block_maxid)   
                blocks[i].append(idx)        
            blocks[i] = tuple(blocks[i])
                    
    return flag, blocks, block_maxid
     


def rank_dict_insert(rank_dict, block_set, rank):
    
    size = len(block_set)
    if not (size in rank_dict):
        rank_dict[size] = {rank:[block_set]}
    elif not (rank in rank_dict[size]):
        rank_dict[size][rank] = [block_set]
    else:        
        rank_dict[size][rank].append(block_set)    
    
    return rank_dict[size][rank]

  
  
def acyclic_detect(block_order):
    
    num = []
    for i in block_order.values():
        num.extend(i)
    in_edge = {i:num.count(i) for i in block_order.keys()}   

    idx = 0
    while idx >= 0:
        idx = -1
        for key, value in in_edge.items():
            if value == 0:
                idx = key         
                break
        if idx >= 0:   
            for j in block_order[idx]:
                in_edge[j] -= 1
            del in_edge[idx]
    
   
    if len(in_edge) > 0:                        
        return False
    else:
        return True
            


def block_order_acyclic(cluster_to_block, block_num):    
   
    block_order = {i:[] for i in range(block_num)}  
    for i in cluster_to_block.values():
        if i[0] != -1 and len(i) > 1:
            bi = i[0][0]
            block_order[bi].extend(i[1:])
    
    for i in range(block_num): 
        block_order[i] = list(set(block_order[i]))
    a = acyclic_detect(block_order)
    
    return a, block_order 

            

def cluster_and_block(block_all, cluster_comp):
 
    len1 = len(cluster_comp)
    cluster_to_block = {i:[-1] for i in range(len1)}
    
    for i in range(len1):
        for j in range(len(block_all)):            
            rem = cluster_comp[i].intersection(block_all[j])
            if len(rem) > 0:  
                if len(rem) < len(block_all[j]):
                    if cluster_to_block[i][0] != -1:
                        return False, [j, rem]
                    else:                         
                        cluster_to_block[i][0] = [j, rem]                        
                else:
                    cluster_to_block[i].append(j)               
 
    return True, cluster_to_block
 
 

def block_set_all(block_set, taxa_comp):        
    
    block_all = []
    a = []    
    for i in block_set: 
        a.extend(list(i))
    set1 = set(taxa_comp).difference(set(a))  
    if len(set1) > 0:              
       block_all = [frozenset(set1)]    
    for i in block_set:
        block_all.append(i)
        
    return block_all




def block_making(incomp, cluster_comp, taxa_comp):         
    
    block_idx = {incomp[0][0]:1, incomp[0][1]:2, incomp[0][2]:3}
    idx_block = {1:[incomp[0][0],1], 2:[incomp[0][1],1], 3: [incomp[0][2],1]}
    block_maxid = 3
    block_set = (1,)
    rank_dict = dict()    
    rank_dict[1] = {0:[[2], [3]]}
    cluster_to_block = dict()
    block_order = dict()
    block_all = []
    min_size = 1   
    max_rank = 0 
    flag = 0
    len1 = len(incomp)
    
    while max_rank < len1:    
        if max_rank == len1 - 1:   
            block_set = [idx_block[i][0] for i in block_set]
            block_all = block_set_all(block_set, taxa_comp)  
            tf, cluster_to_block = cluster_and_block(block_all, cluster_comp)
            
            if tf == False:           
                flag = 1
            else:
                acyclic, block_order = block_order_acyclic(cluster_to_block, len(block_all))       
                if acyclic == True:
                    break
                else:
                    flag = 1
                    
        else:               
            tf, blocks, block_maxid = block_resolve(incomp[max_rank + 1], block_set, idx_block, block_idx, block_maxid)             
            if tf == True:
                max_rank += 1                
            else:
                for i in block_set:
                    block_dict_del(idx_block, block_idx, i)
                for i in range(0, 3):
                    rank_dict_insert(rank_dict, blocks[i], max_rank + 1)   
                flag = 1
             
        if flag == 1:
            min_size = min(rank_dict)    
            max_rank = max(rank_dict[min_size])  
            block_set = rank_dict[min_size][max_rank].pop() 
            flag = 0
            if len(rank_dict[min_size][max_rank]) == 0:
                del rank_dict[min_size][max_rank]        
            if len(rank_dict[min_size]) == 0:   
                del rank_dict[min_size]                   
    
    return block_all, cluster_to_block, block_order, min_size         



def output_block_set(sample_dir, result_dir, taxa_occur):
    
    sample_dir = sample_dir.strip()
    if not sample_dir.endswith('/'):
        sample_dir += '/'
    result_dir = result_dir.strip()
    if not result_dir.endswith('/'):
        result_dir += '/'
    result_file = result_dir + 'output.txt'    
    
         
    tree_file_list, tree_dict_list, tree_taxa_dict, taxa_to_origin  = read_files(sample_dir)         
    tree_num = len(tree_file_list)
    set1 = {i:dict() for i in range(0, tree_num-1)}
    for i in range(0, tree_num-1): 
        set1[i] = {j:[] for j in range(i+1, tree_num)}        
        
        
    with open(result_file, 'w') as f:  
        result ='Tree index'+ '\t\t' + 'Tree file\n' 
        print(result, file=f)
        for i in range(len(tree_file_list)):        
            result =str(i)+ '\t\t\t\t' + str(tree_file_list[i]) + '\t\t'
            print(result, file=f)    
        result ='\n\nTree pair'+ '\t' + 'common taxa' + '\t\t'+'cluster number' + '\t\t'+ 'reticulate node\n' 
        print(result, file=f)      
        
        
    if taxa_occur <= 2:
        for i1 in range(0, tree_num-1):
            for i2 in range(i1+1, tree_num):
                cluster_in, taxa_num, cluster_num = init_file_cluster(tree_dict_list, tree_taxa_dict, [i1,i2], taxa_occur)
                taxa_super_list, taxa_comp_list, cluster_comp_list, state_dict = component_summary(cluster_in)   
        
                reti_num = 0
                reti_each = []
                print('\n For two trees', f"({i1}, {i2}),", 'the block_sets are')
                for i in range(1,len(taxa_comp_list)):                 
                    block_set, cluster_to_block, block_order, num = block_making(state_dict[i], cluster_comp_list[i], taxa_comp_list[i])
                    for j in range(len(block_set)):
                        a = []
                        block_set[j] = list(block_set[j])
                        for k in block_set[j]:
                            a.extend(taxa_super_list[i][k])
                        block_set[j] = set([taxa_to_origin[k] for k in a])                    
                           
                    set1[i1][i2].append(block_set)
                    reti_each.append(num)
                    reti_num += num
                    print('\n',block_set)
                print('\n The total number of reticulate nodes is', reti_num, '\n')
        
                result = '(' + str(i1) + ', ' + str(i2)+')' +'\t\t\t'  
                result +=  str(taxa_num) + '\t\t\t\t'+str(cluster_num) + '\t\t\t\t\t'+ str(reti_num) + '\t\t\t\t'+ str(reti_each) 
                with open(result_file, 'a+') as f:  
                    print(result, file=f)

        with open(result_file, 'a+') as f:  
                result = '\n\nThe block set for each biconnected compoment in a network is listed as follows:\n'
                result += '(All taxa in a block set can be included by an element in the block set for an upper component)\n'
                print(result, file=f)
        for i1 in range(0, tree_num-1):
            for i2 in range(i1+1, tree_num):
                with open(result_file, 'a+') as f:  
                    result = 'Tree (' + str(i1) + ', ' + str(i2)+'):'
                    print(result, file=f)
                    result =''
                    for i in set1[i1][i2]:
                        result += str(i)+'\n'
                    print(result, file=f)


def main():
    # block_tree.py is under the directory 'D:/my_dir/HGT_program/'
    sample_dir = 'D:/my_dir/HGT_program/hgt_sample/'    # The directory contains newick tree files
    result_dir = 'D:/my_dir/HGT_program/hgt_output/'    # The directory contains the running result output.txt
    taxa_occur  = 2   # Each taxa occurs in at least taxa_occur trees. 
    output_block_set(sample_dir, result_dir, taxa_occur)
    

main()