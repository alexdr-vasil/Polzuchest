import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# Функции
def exponential(t, A, alpha, E0):
    return E0 + A * (1 - np.exp(alpha * t))


def polinomic(t, a, b, c):
    return a + b * np.power(t, c)


def exponential_2(t, A1):
    tau = 461
    global Parameters_exp, alpha
    return exponential(t, *Parameters_exp) + A1 * (1 - np.exp(alpha * (t - tau)))


def polinomic_2(t, p, q, n):
    tau = 461
    global Parameters_pol
    return polinomic(t, *Parameters_pol) + polinomic((t - tau), p, q, n)


def relaxation(t, d, e, f):
    tau = 783
    global Parameters_pol2
    return polinomic_2(t, *Parameters_pol2) - polinomic((t - tau), d, e, f)


def relaxation_e(t, A2):
    tau = 783
    global Parameters_exp2
    return exponential_2(t, *Parameters_exp2) - A2 * (1 - np.exp(alpha * (t - tau)))


# Ввод данных из файлов
with open("time.txt", 'r') as f:
    t_graph = []
    for i in f.readlines():
        t_graph.append(float(i))

with open("E.txt", 'r') as f:
    e_graph = []
    for i in f.readlines():
        e_graph.append(float(i))

with open("time1.txt", 'r') as f:
    t1 = []
    for i in f.readlines():
        t1.append(float(i))

with open("E1.txt", 'r') as f:
    e1 = []
    for i in f.readlines():
        e1.append(float(i))

with open("time2.txt", 'r') as f:
    t2 = []
    for i in f.readlines():
        t2.append(float(i))

with open("E2.txt", 'r') as f:
    e2 = []
    for i in f.readlines():
        e2.append(float(i))

with open("rel_time.txt", 'r') as f:
    t3 = []
    for i in f.readlines():
        t3.append(float(i))

with open("rel_e.txt", 'r') as f:
    e3 = []
    for i in f.readlines():
        e3.append(float(i))

T1 = np.array(t1)
E1 = np.array(e1)

T2 = np.array(t2)
E2 = np.array(e2)

T3 = np.array(t3)
E3 = np.array(e3)

# Участок для одной сигма
Parameters_exp, Delta_exp = curve_fit(exponential, T1, E1, bounds=([0.0, -0.5, 3.5], [3.5, 0.0, 4.5]))
Parameters_pol, Delta_pol = curve_fit(polinomic, T1, E1)
print()
print("Параметры экспоненциальной ф-ции на 1 участке: " + str(Parameters_exp))
print("Параметры полиноминальной ф-ции  на 1 участке:  " + str(Parameters_pol))

A, alpha, E0 = [0.8676114, -0.01219247, 4.03983289]
a, b, c = [-3.72277478, 7.39744444, 0.02627248]
t_1 = 461

# Участок для двух сигма
Parameters_exp2, Delta_exp2 = curve_fit(exponential_2, T2, E2, bounds=([0.0], [5]))
Parameters_pol2, Delta_pol2 = curve_fit(polinomic_2, T2, E2)
print()
print("Параметры экспоненциальной ф-ции на 2 участке: " + str(Parameters_exp2))
print("Параметры полиноминальной  ф-ции на 2 участке: " + str(Parameters_pol2))

# Участок релаксации
Parameters_rel, Delta_rel = curve_fit(relaxation, T3, E3)
Parameters_rel_exp, Delta_rel_exp = curve_fit(relaxation_e, T3, E3)
print()
print("Параметры экспоненциальной ф-ции на участке релаксации: " + str(Parameters_rel_exp))
print("Параметры полиноминальной  ф-ции на участке релаксации: " + str(Parameters_rel))

x1_list = np.linspace(0, 461, 462)
x2_list = np.linspace(461, 782, 782 - 461 + 1)
x3_list = np.linspace(783, 814, 814 - 783 + 1)


# График экспоненциальных кривых
fig = plt.figure(figsize=(12, 8))
plt.plot(t_graph, e_graph, '.', label="Экспериментальные точки")
plt.plot(x1_list, exponential(x1_list, *Parameters_exp), label="Теоретическая кривая на 1 участке")
plt.plot(x2_list, exponential_2(x2_list, *Parameters_exp2), label="Теоретическая кривая на 2 участке")
plt.plot(x3_list, relaxation_e(x3_list, *Parameters_rel_exp), label="Теоретическая кривая релаксации")
axes = plt.gca()
axes.set_ylim([0, 7])
plt.title("Сравнение теоретической кривой №1\n (экспонента)")
plt.grid()
plt.legend(loc='lower right')
plt.savefig("exp.png")
plt.show()


# График полиноминальных кривых
fig1 = plt.figure(figsize=(12, 8))
plt.plot(t_graph, e_graph, '.', label="Экспериментальные точки", )
plt.plot(x1_list, polinomic(x1_list, *Parameters_pol), label="Теоретическая кривая на 1 участке")
plt.plot(x2_list, polinomic_2(x2_list, *Parameters_pol2), label="Теоретическая кривая на 2 участке")
plt.plot(x3_list, relaxation(x3_list, *Parameters_rel), label="Теор. кривая релаксации")
plt.legend(loc='lower right')
axes = plt.gca()
axes.set_ylim([0, 7])
plt.title("Сравнение теоретической кривой №2\n (полином)")
plt.grid()
plt.savefig("pol.png")
plt.show()
