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

def Th_Ph_roh(B, d, hs): 
    #Th = [K], Ph = [гПа], roh = [г/м^3]
    h = hs/1000
    def low_latitudes(h):
        Th = 300.4222 - 6.3533*h + 0.005886*(h**2)
        Ph = 1012.0306 - 109.0338*h + 3.6316*(h**2)
        roh = 19.6542*np.exp(-0.2313*h-0.1122*(h**2) + 0.01351*(h**3) - 0.0005923*(h**4))
        return Th, Ph, roh

    #средние широты
    def middle_latitudes(d, h):
        #лето
        def summer(h):
            Th = 294.9838 - 5.2159*h - 0.07109*(h**2)
            Ph = 1012.8186 - 111.5569*h + 3.8646*(h**2)
            roh = 14.3542*np.exp(-0.4174*h - 0.0229*(h**2) + 0.001007*(h**3))
            return Th, Ph, roh
            #зима
        def winter(h):
            Th = 272.7241 - 3.6217*h - 0.1759*(h**2)
            Ph = 1018.8627 - 124.2954*h + 4.8307*(h**2)
            roh = 3.4742*np.exp(-0.2697*h - 0.03604*(h**2) + 0.0004489*(h**3))
            return Th, Ph, roh
        if d >= 92 and d <= 282:
            Th, Ph, roh = summer(h)
        elif d >= 1 and d <= 91:
            Th, Ph, roh = winter(h)
        elif d >= 283 and d <= 365:
            Th, Ph, roh = winter(h)
        return Th, Ph, roh

    #высокие широты
    def high_latitudes(d, h):
        def summer(h):
            Th = 286.8374 - 4.7805*h - 0.1402*(h**2)
            Ph = 1008.0278 - 113.2494*h + 3.9408*(h**2)
            roh = 8.988*np.exp(-0.3614*h - 0.005402*(h**2) + 0.001995*(h**3))
            return Th, Ph, roh

        def winter(h):
            Th = 257.4345 + 2.3474*h - 1.5479*(h**2) + 0.08473*(h**3)
            Ph = 1010.8828 - 122.2411*h + 4.554*(h**2)
            roh = 1.2319*np.exp(-0.07481*h - 0.0981*(h**2) + 0.00281*(h**3))
            return Th, Ph, roh
        if d >= 92 and d <= 282:
            Th, Ph, roh = summer(h)
        elif d >= 1 and d <= 91:
            Th, Ph, roh = winter(h)
        elif d >= 283 and d <= 365:
            Th, Ph, roh = winter(h)
        return Th, Ph, roh

    if B <= 22:
        Th, Ph, roh = low_latitudes(h)
    elif B > 22 and B <= 45:
        Th, Ph, roh = middle_latitudes(d, h)
    elif B > 45 and B <= 90:
        Th, Ph, roh = high_latitudes(d, h)
    eh = roh*Th/216.7
    return Th, Ph, eh

def proposfere_delay(alpha, B, d, hs):
    Ts, Ps, es = Th_Ph_roh(B, d, hs)
    #Ps давление на пункте, гПа
    #Ts температура на пункте, К
    #es парциальное давление водяного пара на пункте, гПа
    #hd высота сухого слоя, м
    #hw высота влажного слоя, м (равная 11км)
    #alpha угол возвышения спутника над горизонтом, угл.градусы
    #r радиус вектора пункта от центра земли, м
    r = 6371110
    lc = 0.167 - (0.076 + 0.00015*(Ts - 273.15))*np.exp(-0.3*alpha)
    hd = 40136 + 148.72*(Ts - 273.15)
    hw = 11000
    Elev = 1.92/(alpha**2 + 0.6)
    Dd = (1.552*Ps*hd/(Ts*1e5))/((1 - (np.cos(alpha*np.pi/180)/(1+lc*hd/r)))**(1/2)) - Elev
    Dw = (0.07465*es*hw/Ts**2)/((1-(np.cos(alpha*np.pi/180)/(1+lc*hw/r)))**(1/2)) - Elev
    Ttropo = Dd+Dw
    return Ttropo/c

def ion_delay(yday, UT, az, el, lat_0, lon_0, hm, fc, alpha, **kwargs):
    '''defTEC'''
    E = get_elevation(alpha)
    TEC = get_tec(yday, UT, az, el, lat_0, lon_0, **kwargs)
    M = 1/(1-(6371*1e3*np.cos(E)/(6371*1e3 + hm))**2)**(1/2) #функция отображения
    deltay_ion = 1.34/1e7*TEC*M/(fc**2)
    return deltay_ion

def get_signal_CA(ti, CA): #сигнал со спутника модулированный для интегрирования
    fc = 10.23*1e6*154 #частота несущей
    t_bit = 1/(1.023*1e6) #длительность бита
    ibit = int(ti/t_bit) #номер бита
    CA_bit_CA = CA[ibit%len(CA)] #СА сигнал
    signal_mod_CA = np.sin(2*np.pi*fc*ti-np.pi*(CA_bit_CA+1)/2) #модулированный сигнал
    return signal_mod_CA, CA_bit_CA

code_length = 1023 #длина кодовой последовательности
prn=1 #номер спутника
d=10 #день
B=30 #широта, град
hs=550 #высота пункта, м
alpha=90 #угол возвышения спутника

UT=14.5
az=np.pi/2
el=alpha*np.pi/180
lat_0=B*np.pi/180
lon_0=104*np.pi/180
fc=f1 #частота несущей
hm = 450*1e3 #высота максимальной электроной плотности
kargs = {"z_start": 80, "z_end": 1000, "l_step": 10,
                  "ne_0": 2e12, "hmax": 300, "half_thickness": 40}
deltay_ion = ion_delay(d, UT, az, el, lat_0, lon_0, hm, fc, alpha, **kargs)

def get_signal_CA(ti, CA): #сигнал со спутника модулированный для интегрирования
    fc = 10.23*1e6*154 #частота несущей
    t_bit = 1/(1.023*1e6) #длительность бита
    ibit = int(ti/t_bit) #номер бита
    CA_bit_CA = CA[ibit%len(CA)] #СА сигнал
    signal_mod_CA = np.sin(2*np.pi*fc*ti-np.pi*(CA_bit_CA+1)/2) #модулированный сигнал
    return signal_mod_CA, CA_bit_CA

#модулированный сигнал для интегрирования
def get_signal(ti, CA, fc, alpha, B, d, hs): #сигнал со спутника модулированный для интегрирования
    elevation = get_elevation(alpha)
    Ttropo = proposfere_delay(alpha, B, d, hs)
    dop = v*np.cos(elevation)/c
    do = fc*dop
    fc = fc + do #частота несущей
    t_bit = 1/(fCA + fCA*dop) #длительность бита
    ti = ti - Ttropo - deltay_ion
    ibit = int(ti/t_bit) #номер бита
    CA_bit = CA[ibit%len(CA)] #СА сигнал
    signal_mod = np.sin(2*np.pi*fc*ti-np.pi*(CA_bit+1)/2) #модулированный сигнал
    return signal_mod, CA_bit

signal_mod_all_CA = []
digital_signal_S_CA = []
signal_mod_all = []
time =[] #ось времени
digital_signal_S = []
CA = get_code_prn(code_length, prn)
for ti in arange(0/fCA, 11/fCA, 1/fs): #для построения графиков
    time.append(ti)
    signal_mod_CA, CA_bit_CA = get_signal_CA(ti, CA)
    signal_mod_all_CA.append(signal_mod_CA) #сигнал без эффектов
    digital_signal_S_CA.append(CA_bit_CA) #CA без эффектов
    signal_mod, CA_bit = get_signal(ti, CA, fc, alpha, B, d, hs) #
    signal_mod_all.append(signal_mod) #сигнал с эффектами
    digital_signal_S.append(CA_bit) #CA с эффектами

fig = plt.figure()
ax_1 = fig.add_subplot(2, 1, 1) #для сигнала на спутнике
ax_1.plot(time, digital_signal_S_CA)
ax_1.set_xlabel('время (с)')
ax_1.set_ylabel('бинарный код без воздействия', fontsize = 7)
ax_1.set_xlim(4.88/1e6, 4.91/1e6)
ax_2 = fig.add_subplot(2, 1, 2) #для сигнала CA
ax_2.plot(time, digital_signal_S)
ax_2.set_xlabel('время (с)')
ax_2.set_ylabel('бинарный код с воздействием', fontsize = 7)
ax_2.set_xlim(4.88/1e6, 4.91/1e6)
show()