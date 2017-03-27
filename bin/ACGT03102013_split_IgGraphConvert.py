'''
Created on 2012. 8. 7.

@author: Akard3
'''
import os
import sys
import time

sys.setrecursionlimit(3000)

print time.ctime()


ms2db_file_name = sys.argv[1]
out_file = sys.argv[2]
input_length = int(sys.argv[3])

file_count = 1
(file_n,ext_n) = os.path.splitext(out_file)
out_file_name = file_n+"_"+str(file_count)+ext_n

f = open(ms2db_file_name,'r')
s = open(out_file_name,'w')
aa_length = 29
min_exon_len = int(input_length)

prev_seq = {}


'''


f = open('VCF.ms2db','r')
#f = open('tmp_ATCG.ms2db','r')
s = open('VCF_chr19.fa','w')
aa_length = 29
min_exon_len = 30
'''




ForwardCode = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
               "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
               "TAT":"Y", "TAC":"Y", "TAA":"X", "TAG":"X",
               "TGT":"C", "TGC":"C", "TGA":"X", "TGG":"W",
               "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
               "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
               "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
               "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
               "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
               "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
               "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
               "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
               "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
               "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
               "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
               "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
               }
ReverseCode =   {"AAA":"F", "AAG":"F", "AAT":"L", "AAC":"L",
               "AGA":"S", "AGG":"S", "AGT":"S", "AGC":"S",
               "ATA":"Y", "ATG":"Y", "ATT":"X", "ATC":"X",
               "ACA":"C", "ACG":"C", "ACT":"X", "ACC":"W",
               "GAA":"L", "GAG":"L", "GAT":"L", "GAC":"L",
               "GGA":"P", "GGG":"P", "GGT":"P", "GGC":"P",
               "GTA":"H", "GTG":"H", "GTT":"Q", "GTC":"Q",
               "GCA":"R", "GCG":"R", "GCT":"R", "GCC":"R",
               "TAA":"I", "TAG":"I", "TAT":"I", "TAC":"M",
               "TGA":"T", "TGG":"T", "TGT":"T", "TGC":"T",
               "TTA":"N", "TTG":"N", "TTT":"K", "TTC":"K",
               "TCA":"S", "TCG":"S", "TCT":"R", "TCC":"R",
               "CAA":"V", "CAG":"V", "CAT":"V", "CAC":"V",
               "CGA":"A", "CGG":"A", "CGT":"A", "CGC":"A",
               "CTA":"D", "CTG":"D", "CTT":"E", "CTC":"E",
               "CCA":"G", "CCG":"G", "CCT":"G", "CCC":"G",
               }



nu_length = aa_length * 3
count = 0

s2 = open('/home/s3cha/data/SpliceDB/NEW_IG_Graph/temp.txt','w')

def getAA(output_array):
    output_seq = []
    for j in range(3):
        AA_seq = ''
        if ForwardFlag == 0:
            for i in range(len(output_array[j:])/3):
                AA_seq += ForwardCode[output_array[3*i+j:3*(i+1)+j]]
        else:
            for i in range(len(output_array[j:])/3):
                AA_seq += ForwardCode[output_array[3*i+j:3*(i+1)+j]]
        output_seq.append(AA_seq)
    return output_seq

def coord_writer(begin,end,check):
    path = [1,3,5,7,9]
    i = 1
    if i != 0:
        if previous_exon[path[i]][path[i-1]] == 3:
            s.write(':')
        else:
            s.write(';')
    return


def write(output_array, path,coordi_x,coordi_y):
    global s
    global count
    if len(output_array)<30:
        return
    if 'X' in output_array or 'N' in output_array:
        return
    AA_array = getAA(output_array)
    temp_s = []
    temp_e = []
    check_deletion = []
    for i in range(len(path)):
        if path[i] != path[-1]:
            if next_exon[path[i]][path[i+1]] == 3:
                check_deletion.append(path[i])
    check = True
    cor_pr = -1
    for j,i in enumerate(path):
        temp_s.append(coordi_start[i])
        temp_e.append(coordi_end[i])
        if correct_frame[i] != -1 and check:
            check = False
            pr_len = 0
            for index in range(j):
                pr_len += temp_e[index]-temp_s[index]
            cor_pr = (correct_frame[i]-pr_len)%3
                
    if ForwardFlag == 0:
        temp_e[0] = coordi_x
        temp_s[len(temp_s)-1] = coordi_y
        if len(temp_s) == 1:
            temp_s[0] = coordi_y
            temp_e[0] = coordi_x
    else:
        temp_s[0] = coordi_x
        temp_e[len(temp_e)-1] = coordi_y
    if ForwardFlag == 0:
        i = 0
        while i<len(temp_s)-1:
            while temp_s[i] == temp_e[i+1]:
                temp_s[i] = temp_s[i+1]
                del temp_s[i+1]
                del temp_e[i+1]
                if i == len(temp_s)-1:
                    break
            i += 1
    else:
        i = 0
        while i<len(temp_s)-1:
            while temp_e[i] == temp_s[i+1]:
                temp_e[i] = temp_e[i+1]
                del temp_s[i+1]
                del temp_e[i+1]
                if i == len(temp_s)-1:
                    break
            i += 1
    code_len = 0
    s_len = temp_e[0] - temp_s[0]
    e_len = temp_e[len(temp_s)-1] - temp_s[len(temp_s)-1]
    for j in range(len(temp_s)):
        code_len += (temp_e[j] - temp_s[j])
    for j in range(3):
        if not check:
            if j != cor_pr:
#                 print AA_array[j]
                s2.write(AA_array[j]+'\n')
                continue
        s2.write('NOTE: '+AA_array[j]+'\n')
        if prev_seq.has_key(AA_array[j]) or 'X' in AA_array[j]:
            continue
        else:
            prev_seq[AA_array[j]] = 0
        s.write('>Splice@'+geneName+'@'+str(count)+'@'+str(ForwardFlag)+';')
        count = count + 1
        for i in range(len(temp_s)):
            if ForwardFlag == 0:
                if i == 0:
                    if s_len == 1:
                        if j == 0:
                            s.write(str(temp_s[i])+'/'+str(temp_e[i]-j)+';')
                        else:
                            continue
                    elif s_len == 2:
                        if j == 0 or 1:
                            s.write(str(temp_s[i])+'/'+str(temp_e[i]-j)+';')
                        else:
                            continue
                    else:
                        s.write(str(temp_s[i])+'/'+str(temp_e[i]-j)+';')
                elif i == len(temp_s)-2 and e_len == 1 and (code_len-j)%3 == 2:
                    s.write(str(temp_s[i]+1)+'/'+str(temp_e[i])+';')
                    break
                elif i == len(temp_s)-1:
                    s.write(str(temp_s[i]+(code_len-j)%3)+'/'+str(temp_e[i])+';')
                else:
                    if s_len == 1 and j == 2 and i == 1:
                        s.write(str(temp_s[i])+'/'+str(temp_e[i]-1)+';')
                    else:
                        s.write(str(temp_s[i])+'/'+str(temp_e[i])+';')
            else:
                if i == 0:
                    if s_len == 1:
                        if j == 0:
                            s.write(str(temp_s[i])+'/'+str(temp_e[i]-j)+';')
                        else:
                            continue
                    elif s_len == 2:
                        if j == 0 or 1:
                            s.write(str(temp_s[i]+j)+'/'+str(temp_e[i])+';')
                        else:
                            continue
                    else:
                        s.write(str(temp_s[i]+j)+'/'+str(temp_e[i])+';')
                elif i == len(temp_s)-2 and e_len == 1 and (code_len-j)%3 == 2:
                    s.write(str(temp_s[i])+'/'+str(temp_e[i]-1)+';')
                    break
                elif i == len(temp_s)-1:
                    s.write(str(temp_s[i])+'/'+str(temp_e[i]-(code_len-j)%3)+';')
                else:
                    if s_len == 1 and j == 2 and i == 1:
                        s.write(str(temp_s[i]+1)+'/'+str(temp_e[i])+';')
                    else:
                        s.write(str(temp_s[i])+'/'+str(temp_e[i])+';')
    
        
        for i in check_deletion:
            if ForwardFlag == 0:
                s.write('<deletion>'+ str(coordi_start[i])+' ')
            else:
                s.write('<deletion>'+ str(coordi_end[i])+' ')
        #for i in range(len(path)):
        #    s.write(str(path[i])+' ')
        s.write('\n')
        s.write(AA_array[j] + '\n')
        
        global out_file_name
        global file_count
        if (os.path.getsize(out_file_name)>100*1024*1024):
            global file_n
            global ext_n
            file_count += 1
            out_file_name = file_n+"_"+str(file_count)+ext_n
            s.close()
            s = open(out_file_name,'w')
        


def MPS2(node_index, L_temp, path, temp_array,coordi_x,check_walk):
    global nu_length
    check_walk[len(check_walk)-1] += exon_length[node_index] 
    if node_index in path:
        print 'loop exist in: ' + geneName + ' at path ',
        for i in path:
            print str(i)+' ',
        print node_index
        return path
    output_array = ''
    temp_path = []
    if L_temp - exon_length[node_index] < 0  or node_index > (max_depth+1)*depth_cutoff:
        path.append(node_index)
        for i in range(len(path)):
            if i != len(path)-1:
                output_array += temp_array[path[i]]
                if output_array != '':
                    temp_path.append(path[i])
            else:
                output_array = output_array + temp_array[path[i]][:L_temp]
                temp_path.append(path[i])
        if ForwardFlag == 0:
            coordi_y = coordi_end[node_index] - L_temp
        else:
            coordi_y = coordi_start[node_index] + L_temp
        write(output_array, temp_path,coordi_x,coordi_y)
        del path[len(path)-1]
        return
    elif len(next_exon[node_index]) == 0  or node_index > (max_depth+1)*depth_cutoff:
        path.append(node_index)
        for i in range(len(path)):
            output_array += temp_array[path[i]]
            if output_array != '':
                temp_path.append(path[i])
        if ForwardFlag == 0:
            coordi_y = coordi_start[node_index]
        else:
            coordi_y = coordi_end[node_index]
        write(output_array, temp_path,coordi_x,coordi_y)
        del path[len(path)-1]
        return
    else:
        L_temp = L_temp - exon_length[node_index]
        temp = next_exon[node_index].keys()
        temp.sort()
        for k,j in enumerate(temp):
            if k == 0:
                if check_walk[-1] < min_exon_len and coordi_end[j] < 0:
                    continue
                check_walk.append(check_walk[len(check_walk)-1])
                if next_exon[node_index][j] != 0:
                    check_walk[len(check_walk)-1] = 0
                path.append(node_index)
                MPS2(j, L_temp, path, temp_array,coordi_x,check_walk)
                del check_walk[len(check_walk)-1]
                del path[len(path)-1]
            else:
                check_walk.append(check_walk[len(check_walk)-1])
                if next_exon[node_index][j] != 0 and check_walk[len(check_walk)-1] < min_exon_len:
                    del check_walk[len(check_walk)-1]
                    continue
                elif next_exon[node_index][j] != 0:
                    check_walk[len(check_walk)-1] = 0     
                path.append(node_index)
                sum = nu_length
                temp_array = exon_array[:]
                for i in range(len(path)):
                    sum = sum - exon_length[path[len(path)-i-1]]
                    if sum <= 0 and sum + exon_length[path[len(path)-i-1]] >0:
                        temp_array[path[len(path)-i-1]] = temp_array[path[len(path)-i-1]][-sum:]    
                        if ForwardFlag == 0:
                            coordi_x = coordi_start[path[len(path)-i-1]] +(sum + exon_length[path[len(path)-i-1]]) 
                        else:
                            coordi_x = coordi_end[path[len(path)-i-1]] -(sum + exon_length[path[len(path)-i-1]])
                    elif sum <= 0:
                        temp_array[path[len(path)-i-1]] = ''
                MPS2(j, L_temp, path, temp_array,coordi_x,check_walk)
                del check_walk[len(check_walk)-1]
                del path[len(path)-1]



def MPS(node_index,path,temp_array,coordi_x,check_walk):
    global nu_length
    check_walk[len(check_walk)-1] += exon_length[node_index]
    output_array = ''
    temp_path = []
    if incoming_indicator[node_index] == 2:
        MPS2(node_index, nu_length, path, temp_array, coordi_x,check_walk)
    elif len(next_exon[node_index]) == 0 or node_index > (max_depth+1)*depth_cutoff:
        incoming_indicator[node_index] = 2
        path.append(node_index)
        for i in range(len(path)):
            output_array += temp_array[path[i]]
            if output_array != '':
                temp_path.append(path[i])
        if ForwardFlag == 0:
            coordi_y = coordi_start[node_index]
        else:
            coordi_y = coordi_end[node_index]
        write(output_array, temp_path,coordi_x,coordi_y)
        del path[len(path)-1]
        return
    else:
        incoming_indicator[node_index] = 2
        temp = next_exon[node_index].keys()
        temp.sort()
        for k,j in enumerate(temp):
            if k == 0:
                if check_walk[-1] < min_exon_len and coordi_end[j] < 0:
                    continue
                check_walk.append(check_walk[len(check_walk)-1])
                if next_exon[node_index][j] != 0:
                    check_walk[len(check_walk)-1] = 0
                path.append(node_index)
                MPS(j,path,temp_array,coordi_x,check_walk)
                del check_walk[len(check_walk)-1]
                del path[len(path)-1]
            else:
                check_walk.append(check_walk[len(check_walk)-1])
                if next_exon[node_index][j] != 0 and check_walk[len(check_walk)-1] < min_exon_len:
                    #print geneName,path,node_index,check_walk
                    del check_walk[len(check_walk)-1]
                    continue
                elif next_exon[node_index][j] != 0:
                    check_walk[len(check_walk)-1] = 0                    
                path.append(node_index)
                sum = nu_length
                temp_array = exon_array[:]
                for i in range(len(path)):
                    sum = sum - exon_length[path[len(path)-i-1]]
                    if sum <= 0 and sum + exon_length[path[len(path)-i-1]]>0:
                        temp_array[path[len(path)-i-1]] = temp_array[path[len(path)-i-1]][-sum:]
                        if ForwardFlag == 0:
                            coordi_x = coordi_start[path[len(path)-i-1]] + sum + exon_length[path[len(path)-i-1]]
                        else:
                            coordi_x = coordi_end[path[len(path)-i-1]] - sum - exon_length[path[len(path)-i-1]]
                    elif sum <= 0:
                        temp_array[path[len(path)-i-1]] = ''
                MPS(j,path,temp_array,coordi_x,check_walk)
                del check_walk[len(check_walk)-1]
                del path[len(path)-1]


for line in f:
    if line.find('<Database') > -1:
        continue
    if line.find('<Gene Name=') > -1:
        #print line ##########################################
        temp = line.split('"')
        geneName = temp[1]
        ExonCount = temp[3]
        ForwardFlag = int(temp[7])
        exon_length   = []
        exon_array    = []
        next_exon     = []
        previous_exon = []
        coordi_start  = []
        coordi_end    = []
        source_exon   = []
        correct_frame = []
        incoming_indicator = []
        for i in range(int(ExonCount)):
            exon_length.append(0)
            exon_array.append('')
            next_exon.append({})
            previous_exon.append({})
            incoming_indicator.append(1)
            coordi_start.append(0)
            coordi_end.append(0)
            correct_frame.append(-1)
        continue
    if line.find('<Exon Index') > -1:
        temp = line.split('"')
        current_exon = int(temp[1])
        coordi_start[current_exon] = int(temp[3])
        coordi_end[current_exon]   = int(temp[5])
        if len(temp)>9:
            correct_frame[current_exon] = int(temp[9])
        continue
    if line.find('<ExonSequence') > -1:
        temp = line.split('"')
        exon_length[current_exon] =int(temp[1])
        temp1 = temp[2].split('>')
        temp2 = temp1[1].split('<')
        exon_array[current_exon] = temp2[0]
        if coordi_end[current_exon] < 0:
            coordi_start[current_exon] = coordi_end[current_exon] - exon_length[current_exon]
        continue
    if line.find('<ExtendsExon Index=') > -1:
        temp = line.split('"')
        previous_exon[current_exon][int(temp[1])] = 0
        continue
    if line.find('<LinkFrom Index=') > -1:
        temp = line.split('"')
        previous_exon[current_exon][int(temp[1])] = 1
        continue
    if line.find('<LinkFromInsertion Index=') > -1:
        temp = line.split('"')
        previous_exon[current_exon][int(temp[1])] = 2
        continue
    if line.find('<LinkFromDeletion Index=') > -1:
        temp = line.split('"')
        previous_exon[current_exon][int(temp[1])] = 3
        continue
    if line.find('<LinkFromMutation Index=') > -1:
        temp = line.split('"')
        previous_exon[current_exon][int(temp[1])] = 4
        continue
    if line.find('</Gene>') ==-1:
        continue
    if line.find('</Gene>') > -1:
        error = 0
        for i in range(len(previous_exon)):
            temp = previous_exon[i].keys()
            for k in temp:
                if k >= len(previous_exon):
                    error = 1
                    print 'Unfeasible link existed in: ' + str(current_exon)
                    break
            if error == 1:
                break
            if previous_exon[i] == {}:
                source_exon.append(i)
            else:
                for j in previous_exon[i].keys(): 
                    #temp = previous_exon[i].keys()
                    next_exon[j][i] = previous_exon[i][j]
        if error == 1:
            continue
        depth_cutoff = 2000
        for max_depth in range(int(ExonCount)/depth_cutoff+1):
            if max_depth != 0:
                for j in range(int(ExonCount)):
                    incoming_indicator[j] = 1
            for depth in range(min(int(ExonCount)-depth_cutoff*max_depth,depth_cutoff)):
                i = depth + depth_cutoff*max_depth
                if incoming_indicator[i] != 1:
                    continue
                if ForwardFlag == 0:
                    coordi_x = coordi_end[i]
                else:
                    coordi_x = coordi_start[i]
                global check_walk
                check_walk = [0]
                MPS(i,[],exon_array,coordi_x,check_walk)
        '''
              
        for i in source_exon:
            if ForwardFlag == 0:
                coordi_x = coordi_end[i]
            else:
                coordi_x = coordi_start[i]
            global check_walk
            check_walk = [0]
            MPS(i,[],exon_array,coordi_x,check_walk)
        '''
        continue
        

print time.ctime()
print 'end'
