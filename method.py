import keras
'''
To generate 3-channel PPI data
'''
import numpy as np


def get_CNN_data_3_channel(protein_A_seq,protein_B_seq):
    def get_channel_1_sequence(seq_list):
        def trans_amio_acid_value(char):
            amio_acid_value = 0
            if char == 'A':
                amio_acid_value = 1 * 12 / 255
            elif char == 'C':
                amio_acid_value = 2 * 12 / 255
            elif char == 'D':
                amio_acid_value = 3 * 12 / 255
            elif char == 'E':
                amio_acid_value = 4 * 12 / 255
            elif char == 'F':
                amio_acid_value = 5 * 12 / 255
                # ---------------------------------
            elif char == 'G':
                amio_acid_value = 6 * 12 / 255
            elif char == 'H':
                amio_acid_value = 7 * 12 / 255
            elif char == 'I':
                amio_acid_value = 8 * 12 / 255
            elif char == 'K':
                amio_acid_value = 9 * 12 / 255
            elif char == 'L':
                amio_acid_value = 10 * 12 / 255
                # ---------------------------------
            elif char == 'M':
                amio_acid_value = 11 * 12 / 255
            elif char == 'N':
                amio_acid_value = 12 * 12 / 255
            elif char == 'P':
                amio_acid_value = 13 * 12 / 255
            elif char == 'Q':
                amio_acid_value = 14 * 12 / 255
            elif char == 'R':
                amio_acid_value = 15 * 12 / 255
                # ---------------------------------
            elif char == 'S':
                amio_acid_value = 16 * 12 / 255
            elif char == 'T':
                amio_acid_value = 17 * 12 / 255
            elif char == 'V':
                amio_acid_value = 18 * 12 / 255
            elif char == 'W':
                amio_acid_value = 19 * 12 / 255
            elif char == 'Y':
                amio_acid_value = 20 * 12 / 255
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
        C_num = 0
        D_num = 0
        E_num = 0
        F_num = 0

        G_num = 0
        H_num = 0
        I_num = 0
        K_num = 0
        L_num = 0

        M_num = 0
        N_num = 0
        P_num = 0
        Q_num = 0
        R_num = 0

        S_num = 0
        T_num = 0
        V_num = 0
        W_num = 0
        Y_num = 0

        for sta in seq_list:
            if sta == 'A':
                A_num = A_num + 1
            elif sta == 'C':
                C_num = C_num + 1
            elif sta == 'D':
                D_num = D_num + 1
            elif sta == 'E':
                E_num = E_num + 1
            elif sta == 'F':
                F_num = F_num + 1
            # ---------------------------------
            elif sta == 'G':
                G_num = G_num + 1
            elif sta == 'H':
                H_num = H_num + 1
            elif sta == 'I':
                I_num = I_num + 1
            elif sta == 'K':
                K_num = K_num + 1
            elif sta == 'L':
                L_num = L_num + 1
            # ---------------------------------
            elif sta == 'M':
                M_num = M_num + 1
            elif sta == 'N':
                N_num = N_num + 1
            elif sta == 'P':
                P_num = P_num + 1
            elif sta == 'Q':
                Q_num = Q_num + 1
            elif sta == 'R':
                R_num = R_num + 1
            # ---------------------------------
            elif sta == 'S':
                S_num = S_num + 1
            elif sta == 'T':
                T_num = T_num + 1
            elif sta == 'V':
                V_num = V_num + 1
            elif sta == 'W':
                W_num = W_num + 1
            elif sta == 'Y':
                Y_num = Y_num + 1
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
            elif char == 'C':
                amio_acid_value = 2
            elif char == 'D':
                amio_acid_value = 3
            elif char == 'E':
                amio_acid_value = 4
            elif char == 'F':
                amio_acid_value = 5
                # ---------------------------------
            elif char == 'G':
                amio_acid_value = 6
            elif char == 'H':
                amio_acid_value = 7
            elif char == 'I':
                amio_acid_value = 8
            elif char == 'K':
                amio_acid_value = 9
            elif char == 'L':
                amio_acid_value = 10
                # ---------------------------------
            elif char == 'M':
                amio_acid_value = 11
            elif char == 'N':
                amio_acid_value = 12
            elif char == 'P':
                amio_acid_value = 13
            elif char == 'Q':
                amio_acid_value = 14
            elif char == 'R':
                amio_acid_value = 15
                # ---------------------------------
            elif char == 'S':
                amio_acid_value = 16
            elif char == 'T':
                amio_acid_value = 17
            elif char == 'V':
                amio_acid_value = 18
            elif char == 'W':
                amio_acid_value = 19
            elif char == 'Y':
                amio_acid_value = 20
            else:
                amio_acid_value = -1

            return amio_acid_value
        value_list = []
        for i in range(len(seq_list)-1):
            value_list.append((trans_amio_acid_value(seq_list[i])+(trans_amio_acid_value(seq_list[i+1])-1)*20)/400)

        value_list.append(0)
        assert len(seq_list) == len(value_list), "error: wanted:{0} but got:{1}".format(len(seq_list),
                                                                                           len(value_list))

        return value_list

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