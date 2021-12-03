import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# Функции
def exponential(t, A, alpha, E0):
    return E0 + A * (1 - np.exp(alpha * t))


def polinomic(t, a, b, c):
    return a + b * np.power(t, c)


def exponential_2(t, A1, E1):
    tau = 461
    global Parameters_exp, alpha
    return exponential(t, *Parameters_exp) + E1 + A1 * (1 - np.exp(alpha * (t - tau)))


def polinomic_2(t, p, q, n):
    tau = 461
    global Parameters_pol
    return polinomic(t, *Parameters_pol) + polinomic((t - tau), p, q, n)


def relaxation(t):
    tau = 783
    global Parameters_pol2
    return polinomic_2(t, *Parameters_pol2) - 1.2*polinomic((t - tau), *Parameters_pol)


def relaxation_e(t):
    tau = 783
    global Parameters_exp2, A
    return exponential_2(t, *Parameters_exp2) - 1.12* A * (1 - np.exp(alpha * (t - tau))) - 1.2*E0


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
Parameters_exp2, Delta_exp2 = curve_fit(exponential_2, T2, E2, bounds=([0.0, 0.3], [2, 1.0]))
Parameters_pol2, Delta_pol2 = curve_fit(polinomic_2, T2, E2)
print()
print("Параметры экспоненциальной ф-ции на 2 участке: " + str(Parameters_exp2))
print("Параметры полиноминальной  ф-ции на 2 участке: " + str(Parameters_pol2))


x1_list = np.linspace(0, 461, 462)
x2_list = np.linspace(461, 782, 782 - 461 + 1)
x3_list = np.linspace(783, 850, 850 - 783 + 1)


# График экспоненциальных кривых
fig = plt.figure(figsize=(8, 6))
plt.plot(t_graph, e_graph, '.', label="Экспериментальные точки")
plt.plot(x1_list, exponential(x1_list, *Parameters_exp), label="Теоретическая кривая на 1 участке")
plt.plot([0, 0], [0, exponential(0, *Parameters_exp)], linestyle='--', dashes=(5, 5), color = 'orange')
plt.plot(x2_list, exponential_2(x2_list, *Parameters_exp2), label="Теоретическая кривая на 2 участке")
plt.plot([461, 461], [5, exponential_2(461, *Parameters_exp2)], linestyle='--', dashes=(5, 1), color = 'green')
plt.plot(x3_list, relaxation_e(x3_list), label="Теоретическая кривая релаксации")
plt.plot([783, 783], [6.25, relaxation_e(783)] , linestyle='--', dashes=(5, 5), color = 'red' )
axes = plt.gca()
axes.set_ylim([0, 7])
plt.xlabel("t, с")
plt.ylabel("E, мм")
plt.title("Сравнение теоретической кривой №1 с экспериментом\n (экспонента)")
plt.grid()
plt.legend(loc='lower center')
plt.savefig("exp.png")
plt.show()


# График полиноминальных кривых
fig1 = plt.figure(figsize=(8, 6))
plt.plot(t_graph, e_graph, '.', label="Экспериментальные точки", )
plt.plot(x1_list, polinomic(x1_list, *Parameters_pol), label="Теор. кривая на 1 участке")
plt.plot(x2_list, polinomic_2(x2_list, *Parameters_pol2), label="Теор. кривая на 2 участке")
plt.plot([461, 461], [5, polinomic_2(461, *Parameters_pol2)], color = 'green')
y_rel = relaxation(x3_list)
y_rel[0] = 6.25
plt.plot(x3_list, y_rel, label="Теор. кривая релаксации")
plt.legend(loc='lower center')
axes = plt.gca()
axes.set_ylim([0, 7])
plt.xlabel("t, с")
plt.ylabel("E, мм")
plt.title("Сравнение теоретической кривой №2 с экспериментом \n (полином)")
plt.grid()
plt.savefig("pol.png")
plt.show()



# Сравнение на участке 1
fig3 = plt.figure(figsize=(8, 6))
plt.plot(t_graph[0:43], e_graph[0:43], '.', label="Экспериментальные точки")
plt.plot(x1_list, exponential(x1_list, *Parameters_exp), label="Экспонента")
plt.plot(x1_list, polinomic(x1_list, *Parameters_pol), label="Полином")
axes = plt.gca()
axes.set_ylim([4, 5])
plt.xlabel("t, с")
plt.ylabel("E, мм")
plt.title("Сравнение двух методов аппроксимации \n на участке 1")
plt.grid()
plt.legend(loc='lower center')
plt.savefig("exp1.png")
plt.show()



# Сравнение на участке 2
fig4 = plt.figure(figsize=(8, 6))
plt.plot(t_graph[45:82], e_graph[45:82], '.', label="Экспериментальные точки")
plt.plot(x2_list, exponential_2(x2_list, *Parameters_exp2), label="Экспонента")
plt.plot(x2_list, polinomic_2(x2_list, *Parameters_pol2), label="Полином")
axes = plt.gca()
axes.set_ylim([5.3, 6.3])
plt.xlabel("t, с")
plt.ylabel("E, мм")
plt.title("Сравнение двух методов аппроксимации \n на участке 2")
plt.grid()
plt.legend(loc='lower center')
plt.savefig("exp2.png")
plt.show()


# График экспоненциальных кривых
fig5 = plt.figure(figsize=(4, 10))
plt.plot(t_graph[83::], e_graph[83::], '.', label="Экспериментальные точки")
plt.plot(x3_list, relaxation_e(x3_list), label="Экспонента")
plt.plot(x3_list, y_rel, label="Полином")
axes = plt.gca()
axes.set_ylim([0, 7])
plt.xlabel("t, с")
plt.ylabel("E, мм")
plt.title("Сравнение \n на участке релаксации")
plt.grid()
plt.legend(loc='upper right')
plt.savefig("exp3.png")
plt.show()
