from scipy import integrate
from pylab import*
import numpy as np
import matplotlib.pyplot as plt

                                        #ПОСТОЯННЫЕ
n = 4.4647/1e10
f_0 = 10.23*1e6*(1+n) #номинальная частота c учетом релятивистского влияния
fCA = f_0/10 #частота СА
f1 = f_0*154 #частота L1
f2 = f_0*120 #частота L2
fs = 1.0*1e11 #
v = 3888.88888 #скорость спутника м/с
c = 3*1e8 #скорость света

#ГЕНЕРАЦИЯ СИГНАЛА
#генерация СА кода
def get_code(code_length, bits): 
    a0 = []
    b0 = []
    a = [1] * 10
    b = [1] * 10
    CA = np.ones(code_length)
    for j in range(code_length):
        a0 = a[9] #выводим
        a1 = (a[2]+a[9])%2 #считаем
        a = a[-1:] + a[0:9] #сдвигаем
        a[0] = a1 #вставляем
        b0 = (b[bits[0]]+b[bits[1]])%2 #считаем и выводим для PRN7 спутника
        b1 = (b[5]+b[7]+b[8]+b[9])%2 #считаем
        b2 = (b[1]+b[2] + b1)%2 #считаем
        b = b[-1:] + b[0:9] #сдвигаем
        b[0] = b2 #вставляем
        CA[j] = (a0+b0)%2
        CA[np.where(CA < 1)] = -1
    return CA

#спутники                                           
def get_code_prn(code_length, prn):
    bits = (1,5), (2,6), (3,7), (4,8), (0,8), (1,9), (0,7), (1,8), (2,9), (1,2), (2,3), (4,5), (5,6), (6,7), (7,8), (8,9), (0,3), (1,4), (2,5), (3,6), (4,7), (5,8), (0,2), (3,5), (4,6), (5,7), (6,8), (7,9), (0,5), (1,6), (2,7), (3,8) # make for prh
    return get_code(code_length, bits[prn-1]) #возвращает СА код

#угол возвышения спутника
def get_elevation(alpha): #prn, ti
    #MAKE ccalculations
    return alpha*np.pi/180 #возвращает угол в радианах



#модулированный сигнал для интегрирования
def get_signal(ti, CA, fc, alpha): #сигнал со спутника модулированный для интегрирования
    elevation = get_elevation(alpha)
    #delay_tro = proposfere_delay(B, d, hs, alpha)
    dop = v*np.cos(elevation)/c
    fc = fc + fc*dop #частота несущей
    t_bit = 1/(fCA + fCA*dop) #длительность бита
    ti = ti #delay_tro
    ibit = int(ti/t_bit) #номер бита
    CA_bit = CA[ibit%len(CA)] #СА сигнал
    signal_mod = np.sin(2*np.pi*fc*ti-np.pi*(CA_bit+1)/2) #модулированный сигнал
    return signal_mod

#ДЕМОДУЛЯЦИЯ
#детектор фазы
def detector_phase(code_length, prn, fc, wch, wlen, step1, alpha):
    CA = get_code_prn(code_length, prn)
    t = []
    ca_bit_phase_t = []
    position_0 = 0
    for i in range(1, code_length+3): #для области интегрирования около примерной смены фазы
        chunck_i_position = position_0 + 1/fCA*i - 1/fc*wch/2 #начало области
        chunck_length = 1/fc*wch #длина бласти                                              
        length_window_all = 1/fc*wlen #длина окон интегрирования                             
        step = (1/fc)* step1 #шаг окон интегрирования                                                         
        for j in range(0, int((chunck_length - length_window_all)/step)): #область окон интегрирования
            s = chunck_i_position + step*j #начало 1 окна
            m = chunck_i_position + step*j + length_window_all/2 #конец 1, начало 2 окон
            f = chunck_i_position + step*j + length_window_all #конец 2 окна
            integ_window_1 = integrate.quad(get_signal, s, m, args=((CA, fc, alpha))) #интегрирование 1 окна от s до m
            integ_window_2 = integrate.quad(get_signal, m, f, args=((CA, fc, alpha))) #интегрирование 2 окна от m до f
            if integ_window_1[0]*integ_window_2[0] > 0: #условие для фазы
                position_0 = f - (1/fCA)*i  # сдвиг начальной позиции (1/fc)
                if integ_window_1[0] > 0 and integ_window_2[0] > 0: #условие для знака бита (1 или -1)
                    ca_bit_phase = 1
                elif integ_window_1[0] < 0 and integ_window_2[0] < 0:
                    ca_bit_phase = -1
                ca_bit_phase_t.append(ca_bit_phase) #создание массива битов (1 и -1)
                t.append(f - (1/fc)/2) #время, где меняется фаза
                break
    return  t, ca_bit_phase_t, CA

#восстановление кода
def restore_code(code_length, prn, fc, wch, wlen, step1, alpha):
    t, ca_bit_phase_t, CA = detector_phase(code_length, prn, fc, wch, wlen, step1, alpha)
    quantity_bit_all = []
    ca_bit = [] #востановленный код

    def get_CA_bit(t, fCA): #функция определения кол-во битов между фазами
        t = [0] + t #добавление начала времени
        for j in range(0, len(t)-1):
            delta_t = t[j + 1] - t[j] #время между фазами
            quantity_bit = round(delta_t/(1/fCA)) #кол-во битов между фазами
            quantity_bit_all.append(quantity_bit) #массив кол-ва битов между фазами
        return quantity_bit_all
    quantity_bit_all = get_CA_bit(t, fCA) #вызываем функцию определения кол-во битов

    for k in range(0, len(ca_bit_phase_t)-1):
        ca_bit += sorted([ca_bit_phase_t[k]]*quantity_bit_all[k])
    ca_bit_1 = ca_bit[-1]*(-1)
    while len(ca_bit) < len(CA):
        ca_bit.append(ca_bit_1)
    return t, ca_bit
#t_i - массив вос-х фаз (во временной оси)
#ca_bit - массив вос-х битов

#битовая ошибка
def error_CA(code_length, prn, fc, wch, wlen, step1, alpha):
    t, ca_bit = restore_code(code_length, prn, fc, wch, wlen, step1, alpha)
    CA = get_code_prn(code_length, prn)
    error_bit = 0
    for i in range(0, len(ca_bit)):
        if CA[i] != ca_bit[i]:
            error_bit +=1
    return error_bit, CA, ca_bit
#error_bit - кол-во бит-х ошибок
#CA - массив эт-х элементов кода
#ca_bit - массив вос-х элементов кода

#ошибки
def get_tested_phase_detector(code_length, prn, fc, wch, wlen, step1, alpha):
    error_bit, CA, ca_bit = error_CA(code_length, prn, fc, wch, wlen, step1, alpha)
    #correct, wrong, t, TCAd, shortage, waste = error_phase(code_length, prn, fc, wch, wlen, step1, alpha)
    return {'correct_code': (CA[-1], CA[-2], CA[-3], CA[-4], CA[-5], len(CA)),
            'restored_code': (ca_bit[-1], ca_bit[-2], ca_bit[-3], ca_bit[-4], ca_bit[-5], len(ca_bit)),
            'error bit': error_bit/code_length,
            #'wrong phase': shortage/len(TCAd),
            #'waste phase': waste/len(TCAd),
            #'correct_phase': len(TCAd),
            #'restored_phase': len(t),
            }

code_length = 1023 #длина кодовой последовательности
prn=9 #номер спутника
fc=f1 #частота несущей
wch=6 #длина области для интегрирования
wlen=3 #длина окон интегрирования
step1=1/5 #шаг окон для интегрирования
alpha=60 #угол возвышения спутника

#phase - кол-во несовпадений сдвиг фаз
#t_i - точки времени сдвигов фаз restored_code
#TCAd - точки времени сдвигов фаз correct_code
#error_bit - кол-во несовпадений бит
#CA - СА код
# #ca_bit - восстановленный СА код
print('alpha', alpha, get_tested_phase_detector(code_length, prn, fc, wch, wlen, step1, alpha))