import keras
'''
To generate 3-channel PPI data
'''
import numpy as np


def get_CNN_data_3_channel(protein_A_seq,protein_B_seq):
    def _2_mer_trans(protein_seq):
        out_seq = []
        for i in range(len(protein_seq)-1):
            if protein_seq[i:i+2] == 'AA':
                out_seq.append('A')
            elif protein_seq[i:i+2] == 'AT':
                out_seq.append('B')
            elif protein_seq[i:i+2] == 'AC':
                out_seq.append('C')
            elif protein_seq[i:i+2] == 'AG':
                out_seq.append('D')

            elif protein_seq[i:i+2] == 'TA':
                out_seq.append('E')
            elif protein_seq[i:i+2] == 'TT':
                out_seq.append('F')
            elif protein_seq[i:i+2] == 'TC':
                out_seq.append('G')
            elif protein_seq[i:i+2] == 'TG':
                out_seq.append('H')

            elif protein_seq[i:i+2] == 'CA':
                out_seq.append('I')
            elif protein_seq[i:i+2] == 'CT':
                out_seq.append('J')
            elif protein_seq[i:i+2] == 'CC':
                out_seq.append('K')
            elif protein_seq[i:i+2] == 'CG':
                out_seq.append('L')

            elif protein_seq[i:i+2] == 'GA':
                out_seq.append('M')
            elif protein_seq[i:i+2] == 'GT':
                out_seq.append('N')
            elif protein_seq[i:i+2] == 'GC':
                out_seq.append('O')
            elif protein_seq[i:i+2] == 'GG':
                out_seq.append('P')

        return out_seq

    def get_channel_1_sequence(seq_list):
        def trans_amio_acid_value(char):
            amio_acid_value = 0
            if char == 'A':
                amio_acid_value = 1*9/150
            elif char == 'B':
                amio_acid_value = 2*9/150
            elif char == 'C':
                amio_acid_value = 3*9/150
            elif char == 'D':
                amio_acid_value = 4*9/150
            elif char == 'E':
                amio_acid_value = 5*9/150
            elif char == 'F':
                amio_acid_value = 6*9/150
                # ---------------------------------
            elif char == 'G':
                amio_acid_value = 7*9/150
            elif char == 'H':
                amio_acid_value = 8*9/150
            elif char == 'I':
                amio_acid_value = 9*9/150
            elif char == 'J':
                amio_acid_value = 10*9/150
            elif char == 'K':
                amio_acid_value = 11*9/150
            elif char == 'L':
                amio_acid_value = 12*9/150
                # ---------------------------------
            elif char == 'M':
                amio_acid_value = 13*9/150
            elif char == 'N':
                amio_acid_value = 14*9/150
            elif char == 'O':
                amio_acid_value = 15*9/150
            elif char == 'P':
                amio_acid_value = 16*9/150
                # ---------------------------------
            else:
                amio_acid_value = -1

            return amio_acid_value
        channel_data = []
        for i in range(len(seq_list)):
            channel_data.append(trans_amio_acid_value(seq_list[i]))
        assert len(seq_list) == len(channel_data),"error: wanted:{0} but got:{1}".format(len(seq_list),len(channel_data))
        return channel_data
    def get_channel_2_statistics(seq_list):
        A_num = 0
        B_num = 0
        C_num = 0
        D_num = 0
        E_num = 0
        F_num = 0
        G_num = 0

        H_num = 0
        I_num = 0
        J_num = 0
        K_num = 0
        L_num = 0
        M_num = 0
        N_num = 0

        O_num = 0
        P_num = 0
        for sta in seq_list:
            if sta == 'A':
                A_num = A_num + 1
            elif sta == 'B':
                B_num = B_num + 1
            elif sta == 'C':
                C_num = C_num + 1
            elif sta == 'D':
                D_num = D_num + 1
            elif sta == 'E':
                E_num = E_num + 1
            elif sta == 'F':
                F_num = F_num + 1
            elif sta == 'G':
                G_num = G_num + 1
                # ---------------------------------
            elif sta == 'H':
                H_num = H_num + 1
            elif sta == 'I':
                I_num = I_num + 1
            elif sta == 'J':
                J_num = J_num + 1
            elif sta == 'K':
                K_num = K_num + 1
            elif sta == 'L':
                L_num = L_num + 1
            elif sta == 'M':
                M_num = M_num + 1
            elif sta == 'N':
                N_num = N_num + 1
                # ---------------------------------
            elif sta == 'O':
                O_num = O_num + 1
            elif sta == 'P':
                P_num = P_num + 1

                # ---------------------------------
        total_data = []
        for char in seq_list:
            if char == 'A':
                total_data.append(A_num/len(seq_list))
            elif char == 'B':
                total_data.append(B_num / len(seq_list))
            elif char == 'C':
                total_data.append(C_num / len(seq_list))
            elif char == 'D':
                total_data.append(D_num / len(seq_list))
            elif char == 'E':
                total_data.append(E_num / len(seq_list))
            elif char == 'F':
                total_data.append(F_num / len(seq_list))
                # ---------------------------------
            elif char == 'G':
                total_data.append(G_num / len(seq_list))
            elif char == 'H':
                total_data.append(H_num / len(seq_list))
            elif char == 'I':
                total_data.append(I_num / len(seq_list))
            elif char == 'J':
                total_data.append(J_num / len(seq_list))
            elif char == 'K':
                total_data.append(K_num / len(seq_list))
            elif char == 'L':
                total_data.append(L_num / len(seq_list))
                # ---------------------------------
            elif char == 'M':
                total_data.append(M_num / len(seq_list))
            elif char == 'N':
                total_data.append(N_num / len(seq_list))
            elif char == 'O':
                total_data.append(O_num / len(seq_list))
            elif char == 'P':
                total_data.append(P_num / len(seq_list))

                # ---------------------------------
        assert len(seq_list) == len(total_data), "通道2计算出现问题：wanted:{0} but got:{1} \n str:{2}".format(len(seq_list),len(total_data),seq_list)
        return total_data
    def get_channel_3_group(seq_list):
        def trans_amio_acid_value(char):
            if char == 'A':
                amio_acid_value = 1
            elif char == 'B':
                amio_acid_value = 2
            elif char == 'C':
                amio_acid_value = 3
            elif char == 'D':
                amio_acid_value = 4
            elif char == 'E':
                amio_acid_value = 5
            elif char == 'F':
                amio_acid_value = 6
            elif char == 'G':
                amio_acid_value = 7
            elif char == 'H':
                amio_acid_value = 8
            elif char == 'I':
                amio_acid_value = 9
            elif char == 'J':
                amio_acid_value = 10
            elif char == 'K':
                amio_acid_value = 11
            elif char == 'L':
                amio_acid_value = 12
            elif char == 'M':
                amio_acid_value = 13
            elif char == 'N':
                amio_acid_value = 14
            elif char == 'O':
                amio_acid_value = 15
            elif char == 'P':
                amio_acid_value = 16
            else:
                amio_acid_value = -1

            return amio_acid_value
        value_list = []
        for i in range(len(seq_list)-1):
            #这里用A+B*16是因为a和b均为1~16，但是AA~AP对应的应该是1~16，即B*16==0，由于b==1，因此需要b-1来满足第一轮B-1==0
            #范围是1~256
            value_list.append((trans_amio_acid_value(seq_list[i])+(trans_amio_acid_value(seq_list[i+1])-1)*16)/256)

        value_list.append(0)
        assert len(seq_list) == len(value_list), "error: wanted:{0} but got:{1}".format(len(seq_list),
                                                                                           len(value_list))

        return value_list

    protein_A_seq = _2_mer_trans(protein_A_seq)
    protein_B_seq = _2_mer_trans(protein_B_seq)
    channal_matrix = np.zeros((3600,3))

    c11 = get_channel_1_sequence(protein_A_seq)
    c12 = get_channel_2_statistics(protein_A_seq)
    c13 = get_channel_3_group(protein_A_seq)

    c21 = get_channel_1_sequence(protein_B_seq)
    c22 = get_channel_2_statistics(protein_B_seq)
    c23 = get_channel_3_group(protein_B_seq)

    for i in range(len(protein_A_seq)):
        channal_matrix[1799 - i][0] = c11[i]
        channal_matrix[1799 - i][1] = c12[i]
        channal_matrix[1799 - i][2] = c13[i]
    for j in range(len(protein_B_seq)):
        channal_matrix[1800 + j][0] = c21[j]
        channal_matrix[1800 + j][1] = c22[j]
        channal_matrix[1800 + j][2] = c23[j]
    return channal_matrix
def get_CFS_data_3_channel(gene_seq):
    '''
    |__1__2__3__4__5__6__7__|  = 8*2+7*1 = 23加边界2 = 25
    :param protein_seq:
    :return:
    '''
    def _2_mer_trans(protein_seq):
        out_seq = []
        for i in range(len(protein_seq)-1):
            if protein_seq[i:i+2] == 'AA':
                out_seq.append('A')
            elif protein_seq[i:i+2] == 'AT':
                out_seq.append('B')
            elif protein_seq[i:i+2] == 'AC':
                out_seq.append('C')
            elif protein_seq[i:i+2] == 'AG':
                out_seq.append('D')

            elif protein_seq[i:i+2] == 'TA':
                out_seq.append('E')
            elif protein_seq[i:i+2] == 'TT':
                out_seq.append('F')
            elif protein_seq[i:i+2] == 'TC':
                out_seq.append('G')
            elif protein_seq[i:i+2] == 'TG':
                out_seq.append('H')

            elif protein_seq[i:i+2] == 'CA':
                out_seq.append('I')
            elif protein_seq[i:i+2] == 'CT':
                out_seq.append('J')
            elif protein_seq[i:i+2] == 'CC':
                out_seq.append('K')
            elif protein_seq[i:i+2] == 'CG':
                out_seq.append('L')

            elif protein_seq[i:i+2] == 'GA':
                out_seq.append('M')
            elif protein_seq[i:i+2] == 'GT':
                out_seq.append('N')
            elif protein_seq[i:i+2] == 'GC':
                out_seq.append('O')
            elif protein_seq[i:i+2] == 'GG':
                out_seq.append('P')

        return out_seq

    def get_channel_1_sequence(seq_list):
        def trans_amio_acid_value(char):
            amio_acid_value = 0
            if char == 'A':
                amio_acid_value = 1*9/150
            elif char == 'B':
                amio_acid_value = 2*9/150
            elif char == 'C':
                amio_acid_value = 3*9/150
            elif char == 'D':
                amio_acid_value = 4*9/150
            elif char == 'E':
                amio_acid_value = 5*9/150
            elif char == 'F':
                amio_acid_value = 6*9/150
                # ---------------------------------
            elif char == 'G':
                amio_acid_value = 7*9/150
            elif char == 'H':
                amio_acid_value = 8*9/150
            elif char == 'I':
                amio_acid_value = 9*9/150
            elif char == 'J':
                amio_acid_value = 10*9/150
            elif char == 'K':
                amio_acid_value = 11*9/150
            elif char == 'L':
                amio_acid_value = 12*9/150
                # ---------------------------------
            elif char == 'M':
                amio_acid_value = 13*9/150
            elif char == 'N':
                amio_acid_value = 14*9/150
            elif char == 'O':
                amio_acid_value = 15*9/150
            elif char == 'P':
                amio_acid_value = 16*9/150
                # ---------------------------------
            else:
                amio_acid_value = -1

            return amio_acid_value
        channel_data = []
        for i in range(len(seq_list)):
            channel_data.append(trans_amio_acid_value(seq_list[i]))
        assert len(seq_list) == len(channel_data),"error: wanted:{0} but got:{1}".format(len(seq_list),len(channel_data))
        return channel_data
    def get_channel_2_statistics(seq_list):
        A_num = 0
        B_num = 0
        C_num = 0
        D_num = 0
        E_num = 0
        F_num = 0
        G_num = 0

        H_num = 0
        I_num = 0
        J_num = 0
        K_num = 0
        L_num = 0
        M_num = 0
        N_num = 0

        O_num = 0
        P_num = 0
        for sta in seq_list:
            if sta == 'A':
                A_num = A_num + 1
            elif sta == 'B':
                B_num = B_num + 1
            elif sta == 'C':
                C_num = C_num + 1
            elif sta == 'D':
                D_num = D_num + 1
            elif sta == 'E':
                E_num = E_num + 1
            elif sta == 'F':
                F_num = F_num + 1
            elif sta == 'G':
                G_num = G_num + 1
                # ---------------------------------
            elif sta == 'H':
                H_num = H_num + 1
            elif sta == 'I':
                I_num = I_num + 1
            elif sta == 'J':
                J_num = J_num + 1
            elif sta == 'K':
                K_num = K_num + 1
            elif sta == 'L':
                L_num = L_num + 1
            elif sta == 'M':
                M_num = M_num + 1
            elif sta == 'N':
                N_num = N_num + 1
                # ---------------------------------
            elif sta == 'O':
                O_num = O_num + 1
            elif sta == 'P':
                P_num = P_num + 1

                # ---------------------------------
        total_data = []
        for char in seq_list:
            if char == 'A':
                total_data.append(A_num/len(seq_list))
            elif char == 'B':
                total_data.append(B_num / len(seq_list))
            elif char == 'C':
                total_data.append(C_num / len(seq_list))
            elif char == 'D':
                total_data.append(D_num / len(seq_list))
            elif char == 'E':
                total_data.append(E_num / len(seq_list))
            elif char == 'F':
                total_data.append(F_num / len(seq_list))
                # ---------------------------------
            elif char == 'G':
                total_data.append(G_num / len(seq_list))
            elif char == 'H':
                total_data.append(H_num / len(seq_list))
            elif char == 'I':
                total_data.append(I_num / len(seq_list))
            elif char == 'J':
                total_data.append(J_num / len(seq_list))
            elif char == 'K':
                total_data.append(K_num / len(seq_list))
            elif char == 'L':
                total_data.append(L_num / len(seq_list))
                # ---------------------------------
            elif char == 'M':
                total_data.append(M_num / len(seq_list))
            elif char == 'N':
                total_data.append(N_num / len(seq_list))
            elif char == 'O':
                total_data.append(O_num / len(seq_list))
            elif char == 'P':
                total_data.append(P_num / len(seq_list))

                # ---------------------------------
        assert len(seq_list) == len(total_data), "通道2计算出现问题：wanted:{0} but got:{1} \n str:{2}".format(len(seq_list),len(total_data),seq_list)
        return total_data
    def get_channel_3_group(seq_list):
        def trans_amio_acid_value(char):
            if char == 'A':
                amio_acid_value = 1
            elif char == 'B':
                amio_acid_value = 2
            elif char == 'C':
                amio_acid_value = 3
            elif char == 'D':
                amio_acid_value = 4
            elif char == 'E':
                amio_acid_value = 5
            elif char == 'F':
                amio_acid_value = 6
            elif char == 'G':
                amio_acid_value = 7
            elif char == 'H':
                amio_acid_value = 8
            elif char == 'I':
                amio_acid_value = 9
            elif char == 'J':
                amio_acid_value = 10
            elif char == 'K':
                amio_acid_value = 11
            elif char == 'L':
                amio_acid_value = 12
            elif char == 'M':
                amio_acid_value = 13
            elif char == 'N':
                amio_acid_value = 14
            elif char == 'O':
                amio_acid_value = 15
            elif char == 'P':
                amio_acid_value = 16
            else:
                amio_acid_value = -1

            return amio_acid_value
        value_list = []
        for i in range(len(seq_list)-1):
            #这里用A+B*16是因为a和b均为1~16，但是AA~AP对应的应该是1~16，即B*16==0，由于b==1，因此需要b-1来满足第一轮B-1==0
            #范围是1~256
            value_list.append((trans_amio_acid_value(seq_list[i])+(trans_amio_acid_value(seq_list[i+1])-1)*16)/256)

        value_list.append(0)
        assert len(seq_list) == len(value_list), "error: wanted:{0} but got:{1}".format(len(seq_list),
                                                                                           len(value_list))

        return value_list

    protein_seq = _2_mer_trans(gene_seq)
    channal_matrix = np.zeros((550,3))  # 550->25*22->7*6->7*5+2
    channal_matrix.reshape((25,22,3))
    channal_matrix = []
    c11 = get_channel_1_sequence(protein_seq)
    c12 = get_channel_2_statistics(protein_seq)
    c13 = get_channel_3_group(protein_seq)

    for i in range(len(protein_seq)):
        x_loc = int(i / 7) * 3 + 3
        y_loc = int(i % 7) * 3 + 3
        channal_matrix[x_loc][y_loc][0] = c11[i]
        channal_matrix[x_loc][y_loc][1] = c12[i]
        channal_matrix[x_loc][y_loc][2] = c13[i]

    return channal_matrix
def __main_plus():
    map_list = np.load("interaction pair.npy")
    # This is an example for this fucntion:
    # map_list:[[1,1,2,2],[1,-1,-1,3],[2,4],[3,5],...]
    # sequence:
    # ['MGANNGKQYG...',
    #  'MAAPASRQVR...'
    #     ...,
    #  'MKRGGRDSDR...',]
    sequence = np.load("sequence.npy")

    all_data = []
    for i in range(len(map_list)):
        s1 = sequence[int(map_list[i][0])].replace('\n','')
        s2 = sequence[int(map_list[i][1])].replace('\n','')
        if (len(s1)<1800) and ((len(s2)<1800)):
            all_data.append(get_CNN_data_3_channel(s1,s2))

    np.save("load data.npy",all_data)

__main_plus()

