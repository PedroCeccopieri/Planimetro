import numpy as np
import time

start = time.time()

r = 500
minDist = 0.01

A = 0

resGauss = 0
resGreen = 0

def main():
    global A

    curva = []

    # Retangulo
    curva.append([50.0,-100.0])
    curva.append([50.0,-50.0])
    curva.append([150.0,-50.0])
    curva.append([150.0,-100.0])

    A = 100 * 50
    print(f"Retangulo:\nArea esperada: {A}")

    printPontos(curva)
    calcular(curva)
    printResultado()

    curva = []

    # Quadrado
    curva.append([50.0,-150.0])
    curva.append([50.0,-50.0])
    curva.append([150.0,-50.0])
    curva.append([150.0,-150.0])

    print("----------------------------------------")

    A = 100**2
    print(f"Quadrado:\nArea esperada: {A}")

    printPontos(curva)
    calcular(curva)
    printResultado()

    curva = []

    # Triangulo
    curva.append([50.0,-150.0])
    curva.append([50.0,-50.0])
    curva.append([150.0,-50.0])

    print("----------------------------------------")

    A = 100**2 / 2
    print(f"Triangulo:\nArea esperada: {A}")

    printPontos(curva)
    calcular(curva)
    printResultado()

    curva = []

    # Trapezio
    curva.append([50.0,-150.0])
    curva.append([50.0,-50.0])
    curva.append([250.0,-50.0])
    curva.append([150.0,-150.0])

    print("----------------------------------------")

    A = (200 + 100) * 100 / 2
    print(f"Trapezio:\nArea esperada: {A}")

    printPontos(curva)
    calcular(curva)
    printResultado()

    curva = []

    n = 6

    # Hexagono 1
    l = 20
    for i in np.arange(0, 2*np.pi - 2*np.pi/n, 2*np.pi/n):
        curva.append([75 + l*np.cos(i),-75 - l*np.sin(i)])

    print("----------------------------------------")
    
    A = 3 * np.sqrt(3) * l**2 / 2
    print(f"Hexagono 1:\nArea esperada: {A}")

    printPontos(curva)
    calcular(curva)
    printResultado()

    curva = []

    l = 100
    for i in np.arange(0, 2*np.pi - 2*np.pi/n, 2*np.pi/n):
        curva.append([75 + l*np.cos(i),-75 - l*np.sin(i)])

    print("----------------------------------------")
    
    A = 3 * np.sqrt(3) * l**2 / 2
    print(f"Hexagono 2:\nArea esperada: {A}")

    printPontos(curva)
    calcular(curva)
    printResultado()

    curva = []


def calcular(curva):
    curva = removeDuplicates(curva)
    curva = normalizeCurve(curva)
    cGreen(curva)
    cGauss(curva)

def printPontos(curva):

    print("\nPontos:\n----------------------------------------")
    for p in curva:
        print(f"({round(p[0],2)},{round(p[1],2)})")
    print("----------------------------------------")

def printResultado():
    print(f"Green: {resGreen}, Erro: {A - resGreen}")
    print(f"Gauss: {resGauss}, Erro: {A - resGauss}")

def getA(theta):
    return (r/2) * np.cos(theta)

def getB(theta):
    return (r/2) * np.sin(theta)

def getX(theta1, theta):
    return getA(theta1) + (r/2) * np.cos(theta)

def getY(theta1, theta):
    return getB(theta1) + (r/2) * np.sin(theta)

def pQuinas(x, y):
	dis = dist(0,0,x,y)
	h = np.sqrt((r/2)**2 - (dis/2)**2)
	return ((x/2) + h*y/dis, (y/2)-h*x/dis)

def anguloVetorEixoX(x, y):
    
    if (y >= 0): return np.arccos(x/(np.sqrt(x*x + y*y)))
    
    if (y < 0): return -np.arccos(x/(np.sqrt(x*x + y*y)))
    
    return 0

def removeDuplicates(arr):
    unique = []
    un = True
    for i in arr:
        for j in unique:
            if (i[0] == j[0] and i[1] == j[1]):
                un = False
                break

        if (un): unique.append(i)
        un = True

    return unique

def normalizeCurve(curva):
    global s
    curva = removeDuplicates(curva)

    newCurva = []
    s = len(curva)

    for i in range(s):

        pa = curva[i]
        pp = curva[(i+1) % s]
        qnt = np.floor(dist(pa[0],pa[1],pp[0],pp[1]) / minDist)
        dx = (pp[0] - pa[0]) / qnt
        dy = (pp[1] - pa[1]) / qnt

        print(f"criando pontos entre: ({round(pa[0],2)},{round(pa[1],2)}) e ({round(pp[0],2)},{round(pp[1],2)})")
        
        for j in np.arange(qnt+1): newCurva.append([pa[0] + dx * j, pa[1] + dy * j])
        
    curva = removeDuplicates(newCurva)
    s = len(curva)
    print("----------------------------------------")
    print(f"Done. A curva esta com {len(curva)} pontos")
    print("----------------------------------------")
    return curva

def calculaPonto(idx, curva, c, d):

    s = len(curva)
    
    p0 = curva[(idx - 1 + s) % s]
    p1 = curva[(idx + s) % s]
    p2 = curva[(idx + 1 + s) % s]
    
    u1 = p2[0] - p1[0]
    v1 = p2[1] - p1[1]
    
    u2 = p1[0] - p0[0]
    v2 = p1[1] - p0[1]
    
    cosTheta1 = -(c*u1+d*v1)/np.sqrt(u1**2 + v1**2)
    cosTheta2 = -(c*u2+d*v2)/np.sqrt(u2**2 + v2**2)

    return (cosTheta1 + cosTheta2)*dist(p1[0],p1[1],p2[0],p2[1])*(r/2)/(2*np.sqrt(c**2 + d**2))

def cGreen(curva):
    global resGreen
    res = 0
    s = len(curva)
    
    for i in range(s):
        p = curva[i]
        theta1 = anguloVetorEixoX(pQuinas(p[0],p[1])[0],pQuinas(p[0],p[1])[1])
        theta2 = anguloVetorEixoX(p[0]-getA(theta1),p[1]-getB(theta1))
        # calculando as coords
        a = getA(theta1)
        b = getB(theta1)
        x = getX(theta1,theta2)
        y = getY(theta1,theta2)
        c = -(y-b)/(r/2)
        d = (x-a)/(r/2)
        res += calculaPonto(i, curva, c, d)

    resGreen = res

def cGauss(curva):
    global resGauss
    res = 0
    
    for i in range(s):
        p = curva[(i - 1 + s) % s]
        q = curva[(i + s) % s]
        res += p[0] * q[1]
        res -= p[1] * q[0]

    resGauss = abs(res)/2

def dist(x1, y1, x2, y2):
    return np.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1))
        
main()
print()
print(f"Tempo de execucao: {time.time() - start} segundos")