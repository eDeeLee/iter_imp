import blosum as bl
import random
import numpy as np
import copy
import sys

def get_amino(seq, index): # indexの位置のアミノ酸を取得
    if (index >= (len(seq))) | (index < 0):
        return '-'
    return seq[index]

def get_inout(seq, index): # 配列の中か外かを判別（必要ないかも？）
    if (index >= (len(seq)-1)) | (index < 1):
        return '-'
    if seq[index] == '-':
        return '-'
    else:
        return '+'
    
def get_score(matrix, amino1, amino2):
    if ((amino1 == '-') | (amino1 == '+')) & ((amino2 == '-') | (amino2 == '+')): # ギャップ-ギャップペナルティの設定
        return 0
    return matrix[amino1][amino2]

def main():
    penalty = -4 # ギャップのペナルティの設定
    matrix = bl.BLOSUM(62, default=penalty) 

    copia = "ILDFHEKLLHPGIQKTTKLFGETYYFPNSQLLIQNIINECSICNLAK"
    MMULV = "LLDFLLHQLTHLSFSKMKALLERSHSPYYMLNRDRTLKNITETCKACAQVN"
    HTLV = "LQLSPAELHSFTHCGQTALTLQGATTTEASNILRSCHACRGGN"
    RSV = "YPLREAKDLHTALHIGPRALSKACNISMQQAREVVQTCPHCNSA"
    SMRV = "IHEATQAHTLHHLNAHTLRLLYKITREQARDIVKACKQCVVAT"
    MMTV = "LESAQESHALHHQNAAALRFQFHITREQAREIVKLCPNCPDWGS"

    seqs = [copia, MMULV, HTLV, RSV, SMRV, MMTV]

    max_steps = 10000 # 最大繰り返し回数
    no_change = 0
    random.seed(int(sys.argv[1]))

    for i in range(max_steps):
        changed = False
        pre_seqs = copy.deepcopy(seqs)
        seqc_num = random.randrange(round(len(seqs)/2)) + 1
        #seqc_num = 1 # 通常のIIMではこちらを使用
        seq_num = random.sample(range(len(seqs)), seqc_num)
        #seq_num = i % len(seqs) # ランダムではなく順番に選ぶ場合（IIM用）
        if i != 0:
            max_len = max([len(seq) for seq in seqs])
            for k in range(len(seqs)): # 全ての配列の長さを揃える（必要ないかも）
                seqs[k] = seqs[k] + '-'*(max_len - len(seqs[k]))

            for k in range(max_len): # SeqCの中で全てがギャップである部分は消す
                space = True
                for j in range(seqc_num):
                    if (seqs[seq_num[j]][max_len - 1 - k] != '-') & (seqs[seq_num[j]][max_len - 1 - k] != '+'):
                        space = False
                        break
                if space:
                    for j in range(seqc_num):
                        seqs[seq_num[j]] = seqs[seq_num[j]][:max_len - 1 - k] + seqs[seq_num[j]][max_len - k:]

            for k in range(max_len): # SeqC以外の中で全てがギャップである部分は消す
                space = True
                for j in range(len(seqs)):
                    if j in seq_num:
                        continue
                    if (seqs[j][max_len - 1 - k] != '-') & (seqs[j][max_len - 1 - k] != '+'):
                        space = False
                        break
                if space:
                    for j in range(len(seqs)):
                        if j in seq_num:
                            continue
                        seqs[j] = seqs[j][:max_len - 1 - k] + seqs[j][max_len - k:]

        seqc = [seqs[num] for num in seq_num]
        seqc_len = max([len(seq) for seq in seqc])
        profile_len = max([len(seq) for seq in seqs])
        DP = np.zeros((seqc_len+1, profile_len+1))
        for j in range(profile_len+1):
            DP[0][j] = penalty * j
        for j in range(seqc_len+1):
            DP[j][0] = penalty * j
        TB = np.zeros((seqc_len+1, profile_len+1)) # トレースバック用

        for j in range(1, seqc_len+1):
            for k in range(1, profile_len+1):
                val1 = 0
                val2 = 0
                val3 = 0
                
                for l in range(len(seqs)):
                    if l in seq_num:
                        continue
                    for m in range(seqc_num):
                        val1 = val1 + get_score(matrix, get_amino(seqc[m], j-1), get_amino(seqs[l], k-1))
                        val2 = val2 + get_score(matrix, get_inout(seqc[m], j), get_amino(seqs[l], k-1))
                        val3 = val3 + get_score(matrix, get_amino(seqc[m], j-1), get_inout(seqs[l], k))
                val1 = (val1 / (len(seqs) - seqc_num)) / seqc_num
                val2 = (val2 / (len(seqs) - seqc_num)) / seqc_num
                val3 = (val3 / (len(seqs) - seqc_num)) / seqc_num

                score1 = DP[j-1][k-1] + val1
                score2 = DP[j][k-1] + val2
                score3 = DP[j-1][k] + val3

                DP[j][k] = max(score1, score2, score3)

                if DP[j][k] == score1: # ギャップは少ない方がいいのでscoreが同じ場合はscore1を優先したい
                    TB[j][k] = 1
                elif DP[j][k] == score2:
                    TB[j][k] = 2
                else:
                    TB[j][k] = 3
        # print(DP[len(seqc)][profile_len])  # 今回のSP score
        s1 = seqc_len
        s2 = profile_len
        while (s1 >= 0) & (s2 >= 0): # トレースバック用のmatrixを見て、適宜ギャップを挿入する
            if TB[s1][s2] == 1:
                s1 = s1 - 1
                s2 = s2 - 1
            elif TB[s1][s2] == 2:
                s2 = s2 - 1
                for l in range(seqc_num):
                    if len(seqs[seq_num[l]]) > s1:
                        if (seqs[seq_num[l]][s1-1] != '-') & (seqs[seq_num[l]][s1] != '-'):
                            seqs[seq_num[l]] = seqs[seq_num[l]][:s1] + '+' + seqs[seq_num[l]][s1:]
                        else:
                            seqs[seq_num[l]] = seqs[seq_num[l]][:s1] + '-' + seqs[seq_num[l]][s1:]
                        #changed = True
            elif TB[s1][s2] == 3:
                s1 = s1 - 1
                
                for l in range(len(seqs)):
                    if l in seq_num:
                        continue
                    if len(seqs[l]) > s2:
                        if (seqs[l][s2-1] != '-') & (seqs[l][s2] != '-'):
                            seqs[l] = seqs[l][:s2] + '+' + seqs[l][s2:]
                        else:
                            seqs[l] = seqs[l][:s2] + '-' + seqs[l][s2:]
                        #changed = True
            elif (s1 == 0) & (s2 > 0):
                s2 = s2 - 1
                for l in range(seqc_num):
                    seqs[seq_num[l]] = '-' + seqs[seq_num[l]]
            elif (s1 > 0) & (s2 == 0):
                s1 = s1 - 1
                for l in range(len(seqs)):
                    if l in seq_num:
                        continue
                    seqs[l] = '-' + seqs[l]
            else:
                break
        
        max_len = max([len(seq) for seq in seqs])
        for k in range(len(seqs)): # 全ての配列の長さを揃える
            seqs[k] = seqs[k] + '-'*(max_len - len(seqs[k]))
    
        for k in range(max_len): # 全ての配列でギャップになっている部分があれば削除
            space = True
            for j in range(len(seqs)):
                if (seqs[j][max_len - 1 - k] != '-') & (seqs[j][max_len - 1 - k] != '+'):
                    space = False
                    break
            if space:
                for j in range(len(seqs)):
                    seqs[j] = seqs[j][:max_len - 1 - k] + seqs[j][max_len - k:]

        if pre_seqs != seqs: # 今回のアラインメントで変化したかどうかを判定
            changed = True

        if changed:
            no_change = 0
        no_change = no_change + 1
        if no_change > 300: # 300回連続で変化がなければ終了
            print('stoped at step %d'%i)
            break

    final_score = 0
    max_len = max([len(seq) for seq in seqs])
    print('random seed is %d'%int(sys.argv[1]))
    for i in range(len(seqs)): # 最終的に得られたマルチプルアラインメントのSP scoreを算出
        print(seqs[i])
        for j in range(i+1, len(seqs)):
            for k in range(max_len):
                final_score = final_score + get_score(matrix, get_amino(seqs[i], k), get_amino(seqs[j], k))
    print('SP score = %d'%final_score)

    
if __name__ == '__main__':
    main()
