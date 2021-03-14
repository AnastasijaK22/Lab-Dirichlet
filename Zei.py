import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt

from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
import newinterface
import sys

def f(x,y):
    return -1 * (np.pi ** 2) * (x**2 + y**2) * np.exp((np.sin(np.pi * x * y)) ** 2) * (2 * np.cos(2 * np.pi * x * y) + (np.sin(2 * np.pi * x * y))**2)

def u(x,y):
    return np.exp((np.sin(np.pi * x * y)) ** 2)

def m1(y):
    return 1

def m2(y):
    return np.exp((np.sin(np.pi * y)) ** 2)

def m3(x):
    return 1

def m4(x):
    return np.exp((np.sin(np.pi * x)) ** 2)

def f_main(x,y):
    return np.sin(np.pi * x * y) ** 2

def m1_main(y):
    return np.sin(np.pi * y)

def m2_main(y):
    return np.sin(np.pi * y)

def m3_main(x):
    return x - x * x

def m4_main(x):
    return x - x * x

def trueMatrix(n,m):
    U = np.zeros((m + 1, n + 1))
    a = 0
    b = 1
    c = 0
    d = 1
    h = (b - a) / n
    k = (d - c) / m
    for j in range(m + 1):
        for i in range(n + 1):
            U[j][i] = u(a + i * h, c + j * k)
    return U

def initialMatrixF(n,m, test):
    F = np.zeros((m + 1, n + 1))
    a = 0
    b = 1
    c = 0
    d = 1
    h = (b - a) / n
    k = (d - c) / m
    if test == 0:
        for j in range(1,m):
            for i in range(1,n):
                F[j][i] = f(a + i * h, c + j * k)
    else:
        for j in range(1,m):
            for i in range(1,n):
                F[j][i] = f_main(a + i * h, c + j * k)
#    print(F)
    return F

def initialMatrixV(n, m, test):
    V = np.zeros((m + 1, n + 1))
    a = 0
    b = 1
    c = 0
    d = 1
    h = (b-a)/n
    k = (d-c)/m
    if test == 0:
        for j in range(m + 1):
            V[j][0] = m1(c + j * k)
        for j in range(m + 1):
            V[j][n] = m2(c + j * k)
        for i in range(n + 1):
            V[0][i] = m3(a + i * h)
        for i in range(n + 1):
            V[m][i] = m4(a + i * h)
    else:
        for j in range(m + 1):
            V[j][0] = m1_main(c + j * k)
        for j in range(m + 1):
            V[j][n] = m2_main(c + j * k)
        for i in range(n + 1):
            V[0][i] = m3_main(a + i * h)
        for i in range(n + 1):
            V[m][i] = m4_main(a + i * h)
    return V

def residual(n, m, V, F):
    res = np.zeros((m + 1, n + 1))
    h_quad = -1 * n * n
    k_quad = -1 * m * m
    A = -2 * (h_quad + k_quad)
    for j in range(1,m):
        for i in range(1,n):
            res[j][i] = A * V[j][i] + h_quad * (V[j][i - 1] + V[j][i + 1])
            res[j][i] = res[j][i] + k_quad * (V[j + 1][i] + V[j - 1][i]) - F[j][i]
    return res

def tsCounter(n, m, V, F, re):
    scal1, scal2 = 0, 0
    h_quad = -1 * n * n
    k_quad = -1 * m * m
    A = -2 * (h_quad + k_quad)
    for j in range(1, m):
        for i in range(1,n):
            up, down, left, right = 1, 1, 1, 1
            if (i == 1):
                left = 0
            if (i == n - 1):
                right = 0
            if (j == 1):
                down = 0
            if (j == m - 1):
                up = 0
            temp = A * re[j][i] + h_quad * (left * re[j][i - 1] + right * re[j][i + 1])
            temp += k_quad * (up * re[j + 1][i] + down * re[j - 1][i])
            scal1 += temp * re[j][i]
            scal2 += temp * temp
    return scal1 / scal2

def methodMR(n, m, SMax, eps, test, S):
    eps_max = 0
    eps_cur = 0
    S[0] = 0
    V = initialMatrixV(n,m, test)
    F = initialMatrixF(n,m, test)
    v_old = 0
    v_new = 0
    flag = False
    while(flag != True):
        eps_max = 0
        re = residual(n,m,V,F)
        ts = tsCounter(n,m,V,F,re)
#        print("re")
#        print(re)
#        print(ts)
#        print("V")
#        print(V)
        temp = 0
        for j in range(1, m):
            for i in range(1, n):
                v_old = V[j][i]
                v_new = V[j][i] - ts * re[j][i]
                temp += re[j][i] * re[j][i]
                eps_cur = np.fabs(v_old - v_new)
                if (eps_cur > eps_max):
                    eps_max = eps_cur
                V[j][i] = v_new
        S[0] = S[0] + 1
        S[2] = np.sqrt(temp)

        if ((eps_max < eps) or (S[0] >= SMax)):
            flag = True
    S[1] = eps_max
    return V

class Example(QMainWindow, newinterface.Ui_MainWindow):
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.pushButton.clicked.connect(self.buttonClicked)

    def buttonClicked(self):

        if (self.lineEdit.text() == ""):
            QMessageBox.question(self, 'Ошибка!', QMessageBox.Ok, QMessageBox.Ok)

        if (self.lineEdit_2.text() == ""):
            QMessageBox.question(self, 'Ошибка!',  QMessageBox.Ok, QMessageBox.Ok)

        if (self.lineEdit_3.text() == ""):
            QMessageBox.question(self, 'Ошибка!',  QMessageBox.Ok, QMessageBox.Ok)

        if (self.lineEdit_4.text() == ""):
            QMessageBox.question(self, 'Ошибка!',  QMessageBox.Ok, QMessageBox.Ok)

        # считываю разбиение по х
        temp = int(self.lineEdit.text())
        N = temp

        # считываю разбиение по у
        temp = int(self.lineEdit_2.text())
        M = temp

        # считываю максимальное число итераций
        temp = int(self.lineEdit_3.text())
        SMax = temp

        temp = float(self.lineEdit_4.text())
        eps = temp

        flag = 0
        if self.checkBox.isChecked() == True:
            flag = 1
#        print(N, M, SMax,eps)

        self.tableWidget.setRowCount(M + 1)
        self.tableWidget.setColumnCount(N + 1)

        a = 0
        b = 1
        c = 0
        d = 1
        h = (b - a) / N
        k = (d - c) / M
#       print(h,k,a,b,c,d)
        S = [0] * 3

        if (flag == 0):

            self.tableWidget_2.setRowCount(M + 1)
            self.tableWidget_2.setColumnCount(N + 1)
            max = 0

            V = methodMR(N, M, SMax, eps, flag, S)
            U = trueMatrix(N, M)
            for i in range(N+1):
                for j in range(M+1):
                    if(max < np.abs(V[j][i] - U[j][i])):
                        max = np.abs(V[j][i] - U[j][i])

            head = []
            for i in range(M, -1, -1):
                temp = "y" + str(i) + " = " + str(a + i * k)
                head.append(temp)

            self.tableWidget.setVerticalHeaderLabels(head)
            self.tableWidget_2.setVerticalHeaderLabels(head)

            head.clear()
            for i in range(0, N+1):
                temp = "x" + str(i) + " = " + str(c + i * h)
                head.append(temp)

            self.tableWidget.setHorizontalHeaderLabels(head)
            self.tableWidget.resizeColumnsToContents()

            self.tableWidget_2.setHorizontalHeaderLabels(head)
            self.tableWidget_2.resizeColumnsToContents()

            for j in range(M, -1, -1):
                for i in range(N + 1):
                    temp = str(V[j][i])
                    newItem = QTableWidgetItem(temp)
                    self.tableWidget.setItem(M - j, i, newItem)
                    temp2 = str(U[j][i])
                    newItem2 = QTableWidgetItem(temp2)
                    self.tableWidget_2.setItem(M - j, i, newItem2)

            resultText = "Для решения тестовой задачи использована сетка с\n"
            resultText += "числом разбиений по х:\n"
            resultText += "n = " + str(N) + "\n"
            resultText += "и числом разбиений по у:\n"
            resultText += "m = " + str(M) + "\n"
            resultText += "применялся метод минимальных невязок\n"
            resultText += "На решение затрачено:\n"
            resultText += "S = " + str(S[0]) + "\n"
            resultText += "Достигнута точность метода: \n"
            resultText += "eps = " + str(S[1]) + "\n"
            resultText += "Задача решена с точностью: \n"
            resultText += "max|vij-uij| = " + str (max)+ "\n"
            resultText += "Евклидова норма невязки:" + str(S[2]) + "\n"

            self.label_9.setText(resultText)
        else:
            V = methodMR(N, M, SMax, eps, flag, S)
            U = trueMatrix(N, M)

            head = []
            for i in range(M, -1, -1):
                temp = "y" + str(i) + " = " + str(a + i * k)
                head.append(temp)

            self.tableWidget.setVerticalHeaderLabels(head)

            head.clear()
            for i in range(0, N + 1):
                temp = "x" + str(i) + " = " + str(c + i * h)
                head.append(temp)

            self.tableWidget.setHorizontalHeaderLabels(head)
            self.tableWidget.resizeColumnsToContents()

            for j in range(M, -1, -1):
                for i in range(N + 1):
                    temp = str(V[j][i])
                    newItem = QTableWidgetItem(temp)
                    self.tableWidget.setItem(M - j, i, newItem)

            resultText = "Для решения основной задачи использована сетка с\n"
            resultText += "числом разбиений по х:\n"
            resultText += "n = " + str(N) + "\n"
            resultText += "и числом разбиений по у:\n"
            resultText += "m = " + str(M) + "\n"
            resultText += "применялся метод минимальных невязок\n"
            resultText += "На решение затрачено:\n"
            resultText += "S = " + str(S[0]) + "\n"
            resultText += "Достигнута точность метода: \n"
            resultText += "eps = " + str(S[1]) + "\n"
            resultText += "Евклидова норма невязки:" + str(S[2]) + "\n"

            self.label_9.setText(resultText)

        # создаём полотно для рисунка
        fig = plt.figure(figsize=(10, 10))

        # создаём рисунок пространства с поверхностью
        ax = fig.add_subplot(1, 1, 1, projection='3d')

        # размечаем границы осей для аргументов
        xval = np.linspace(0, 1, N + 1)
        yval = np.linspace(0, 1, M + 1)

        # создаём массив с xval столбцами и yval строками
        # - в этом массиве будут храниться значения z
        x, y = np.meshgrid(xval, yval)

        # создаём поверхность
        surf = ax.plot_surface(
            # отмечаем аргументы и уравнение поверхности
            x, y, V,
            # шаг прорисовки сетки
            # - чем меньше значение, тем плавнее
            # - будет градиент на поверхности
            # rstride=N + 1,
            # cstride=M + 1,
            cmap=cm.rainbow)

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('v(x,y)')
        temp = "График численного решения"
        ax.set_title(temp)

        fig.canvas.set_window_title('Визуализация решения')
        fig.show()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    form = Example()
    form.show()
    app.exec()