#ввод-вывод с использованием BioPython из текстового файла
#ООП
#тесты
#выводим всю таблицу динамического программирования
#последовательности длиной >10

# 2 скрытых состояния H, L

from Bio import SeqIO
import numpy as np

class Baum_Welch:
    def __init__(self, seq):

        #вероятности
        self.a = {'00': 0.2, '01': 0.8, '10': 0.4, '11': 0.6}
        #0 - H, 1 - L
        self.e = {'A': [0.25, 0.25], 'C': [0.25, 0.25], 'G': [0.25, 0.25], 'T': [0.25, 0.25]}

        self.seq = seq
        self.seq_len = len(seq)
        self.p = 0

    def forward(self):

        print('Просмотр вперед:')
        # инициализация
        self.f_var = 0
        self.f = [[0.0], [0.0]] #первая строка - H, вторая - L, индексы - буквы
        self.summ = 0
        self.letters = {'A': [[-1], [-1]], 'C': [[-1], [-1]], 'G': [[-1], [-1]], 'T': [[-1], [-1]]} #в массивах - индексы, где встречались буквы
        for i in 'ACGT':
            for j in range(self.seq_len):
                self.letters[i][0].append(-1)
                self.letters[i][1].append(-1)
        #print(self.letters)

        # первые элементы отдельно (?), f_0(0) = 1

        self.f_var = self.e[self.seq[0]][0]
        self.summ = self.a['0' + str(0)] + self.a['1' + str(0)]
        self.f[0].append(round(self.f_var * self.summ, 4))

        self.f_var = self.e[self.seq[0]][1]
        self.summ = self.a['0' + str(1)] + self.a['1' + str(1)]
        self.f[1].append(round(self.f_var * self.summ, 4))

        self.letters[self.seq[0]][0][1] = 1
        self.letters[self.seq[0]][1][1] = 1

        # рекурсия
        for i in range(2, self.seq_len + 1):
            for j in range(2):
                self.f_var = self.e[self.seq[i-1]][j]
                self.summ = self.f[0][i-1]*self.a['0'+str(j)] + self.f[1][i-1]*self.a['1'+str(j)]
                self.f[j].append(round(self.f_var * self.summ, 9))
                self.letters[self.seq[i-1]][0][i] = 1
                self.letters[self.seq[i-1]][1][i] = 1
        for k in range(2):
            print('f(', k, '):', self.f[k])

        #for i in 'ACGT':
            #print(self.letters[i])

        self.p = 0
        for i in range(2):
            self.p += self.f[i][self.seq_len - 1]*self.a[str(i)+'0']

        return self.f, self.p, self.letters



    def backward(self):

        print('Просмотр назад:')

        # инициализация
        self.b_var = 0
        self.b = [[self.a['00']], [self.a['10']]]  # первая строка - H, вторая - L, индексы - буквы
        self.summ = 0
        self.letters = {'A': [[-1], [-1]], 'C': [[-1], [-1]], 'G': [[-1], [-1]], 'T': [[-1], [-1]]}  # в массивах - индексы, где встречались буквы
        for i in 'ACGT':
            for j in range(self.seq_len):
                self.letters[i][0].append(-1)
                self.letters[i][1].append(-1)
        #print(self.letters)

        # рекурсия
        for i in range(1, self.seq_len + 1):
            for j in range(2):
                self.b_var = self.e[self.seq[i - 1]][j]
                self.summ = self.b[0][i - 1] * self.a['0' + str(j)] + self.b[1][i - 1] * self.a['1' + str(j)]
                self.b[j].append(round(self.f_var * self.summ, 8))
                self.letters[self.seq[i - 1]][0][i] = 1
                self.letters[self.seq[i - 1]][1][i] = 1

        for k in range(2):
            self.b[k] = list(reversed(self.b[k]))
            print('b(', k, '):', self.b[k])

        for i in 'ACGT':
            self.letters[i][0] = list(reversed(self.letters[i][0]))
            self.letters[i][1] = list(reversed(self.letters[i][1]))
        #print(self.letters)


        self.p = 0
        for i in range(2):
            self.p += self.b[i][0] * self.a['0'+str(i)] * self.e[self.seq[0]][i]

        return self.b, self.p, self.letters


    def new_parametrs(self):


        f, p_forward, letters_f = self.forward()
        b, p_backward, letters_b = self.backward()

        #какую из p нужно использовать?

        #пересчет параметров A_kl
        print('Параметры A_kl:')
        A_kl = [[0, 0], [0, 0]] #[00, 01], [10, 11]
        for k in range(2):
            for l in range(2):
                for i in range(self.seq_len - 1):
                    #print(f[k][i], self.a[str(k)+str(l)], self.e[self.seq[i+1]][l], b[l][i + 1])
                    A_kl[k][l] += round(f[k][i] * self.a[str(k)+str(l)] * self.e[self.seq[i+1]][l] * b[l][i + 1], 9)
                A_kl[k][l] = round(A_kl[k][l]/p_backward, 4)
                print('A(', k, l, ') = ', A_kl[k][l])
        #нормируем
        print('Пересчет + нормировка A_kl:')
        a_kl = [[], []]
        for k in range(2):
            for l in range(2):
                #print(A_kl[k][l], A_kl[k][1-l])
                a_kl[k].append(round(A_kl[k][l]/(A_kl[k][l]+A_kl[k][1-l]), 4))
                print('a(', k, l, ') = ', a_kl[k][l])

        #пересчет параметров E_k(b):
        print('Параметры E_k(b):')
        E_k = {'A': [0, 0], 'C': [0, 0], 'G': [0, 0], 'T': [0, 0]}
        e_k = {'A': [0, 0], 'C': [0, 0], 'G': [0, 0], 'T': [0, 0]}
        summ1 = 0

        for i in 'ACGT':
            #print(i)
            for j in range(1, self.seq_len + 1):
                if letters_f[i][0][j] == 1:
                    E_k[i][0] += f[0][j]*b[0][j]
                    E_k[i][1] += f[1][j]*b[1][j]
            E_k[i][0] = round(E_k[i][0]/p_forward, 6)
            E_k[i][1] = round(E_k[i][1]/p_forward, 6)
            print('E_0(', i, ') = ', E_k[i][0])
            print('E_1(', i, ') = ', E_k[i][1])

        print('Пересчет + нормировка E_k(b):')
        for i in range(2):
            for j in 'ACGT':
                summ1 += E_k[j][i]
            #print(summ1)
            e_k['A'][i] = round(E_k['A'][i] / summ1, 6)
            print('e_', i, '(A)', e_k['A'][i])
            e_k['C'][i] = round(E_k['C'][i] / summ1, 6)
            print('e_', i, '(C)', e_k['C'][i])
            e_k['G'][i] = round(E_k['G'][i] / summ1, 6)
            print('e_', i, '(G)', e_k['G'][i])
            e_k['T'][i] = round(E_k['T'][i] / summ1, 6)
            print('e_', i, '(T)', e_k['T'][i])
            summ1 = 0

if __name__ == '__main__':
    file1 = open("seq1.fasta", "r")
    seq1 = SeqIO.read(file1, "fasta")
    print('Последовательность: ', seq1.seq)

    align = Baum_Welch(seq1.seq)
    align.new_parametrs()
