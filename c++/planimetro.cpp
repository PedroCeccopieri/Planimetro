#include <iostream>
#include <vector>
#include <math.h>
#include <chrono>

using namespace std;

clock_t start = clock();

typedef struct point {
    double x;
    double y;
}point;

void calcular();
void printPontos();
void printResultado();

double getA(double theta);
double getB(double theta);
double getX(double theta1, double theta);
double getY(double theta1, double theta);
double pQuinaX(double x, double y);
double pQuinaY(double x, double y);
double anguloVetorEixoX(double x, double y);
vector<point> removeDuplicates(vector<point> arr);
vector<point> normalizeCurve(vector<point> curva);
double calculaPonto(vector<point> curva, int idx, int s);
void cGreen(vector<point> curva);
void cGauss(vector<point> curva);
double dist(double x1, double y1, double x2, double y2);

double a = 0, b = 0;
double x = 0, y = 0;
double c = 0, d = 0;

double theta1 = 0, theta2 = 0;

double resGauss = 0;
double resGreen = 0;

int r = 500;
double minDist = 0.01;

double A;

vector<point> curva;

int main() {
    
    // Retangulo
    curva.push_back({50.0,-100.0});
    curva.push_back({50.0,-50.0});
    curva.push_back({150.0,-50.0});
    curva.push_back({150.0,-100.0});

    A = 100 * 50;
    printf("Retangulo:\nArea esperada: %f\n",A);

    printPontos();
    calcular();
    printResultado();

    curva.clear();

    // Quadrado
    curva.push_back({50.0,-150.0});
    curva.push_back({50.0,-50.0});
    curva.push_back({150.0,-50.0});
    curva.push_back({150.0,-150.0});

    cout << "----------------------------------------" << endl;

    A = 100 * 100;
    printf("Quadrado:\nArea esperada: %f\n",A);

    printPontos();
    calcular();
    printResultado();

    curva.clear();

    // Triangulo
    curva.push_back({50.0,-150.0});
    curva.push_back({50.0,-50.0});
    curva.push_back({150.0,-50.0});

    cout << "----------------------------------------" << endl;

    A = 100 * 100 / 2;
    printf("Triangulo:\nArea esperada: %f\n",A);

    printPontos();
    calcular();
    printResultado();

    curva.clear();

    // Trapezio
    curva.push_back({50.0,-150.0});
    curva.push_back({50.0,-50.0});
    curva.push_back({250.0,-50.0});
    curva.push_back({150.0,-150.0});

    cout << "----------------------------------------" << endl;

    A = (200 + 100) * 100 / 2;
    printf("Trapezio:\nArea esperada: %f\n",A);

    printPontos();
    calcular();
    printResultado();

    curva.clear();

    int n = 6;
    int l;

    // Hexagono 1
    l = 20;
    for (double i = 0; i < 2*M_PI - 2*M_PI/n; i += 2*M_PI/n) curva.push_back({75 + l*cos(i),-75 - l*sin(i)});

    cout << "----------------------------------------" << endl;

    A = 3 * sqrt(3) * l*l / 2;
    printf("Hexagono 1:\nArea esperada: %f\n",A);

    printPontos();
    calcular();
    printResultado();

    curva.clear();

    // Hexagono 2
    l = 100;
    for (double i = 0; i < 2*M_PI - 2*M_PI/n; i += 2*M_PI/n) curva.push_back({75 + l*cos(i),-75 - l*sin(i)});

    cout << "----------------------------------------" << endl;

    A = 3 * sqrt(3) * l*l / 2;
    printf("Hexagono 2:\nArea esperada: %f\n",A);

    printPontos();
    calcular();
    printResultado();

    curva.clear();

    cout << endl << "Tempo de execucao: " << (double)(clock() - start)/CLOCKS_PER_SEC << " segundos" << endl; 
    
    return 0;
}

void calcular() {
    curva = removeDuplicates(curva);
    curva = normalizeCurve(curva);
    cGreen(curva);
    cGauss(curva);
}

void printPontos() {

    cout << endl << "Pontos:"  << endl << "----------------------------------------" << endl;
    for(auto x:curva) {
        cout << "(" << x.x << "," << x.y << ")\n";
    }
    cout << "----------------------------------------" << endl;
}

void printResultado() {
    printf("Green: %f, Erro: %f\n",resGreen, A - resGreen);
    printf("Gauss: %f, Erro: %f\n",resGauss, A - resGauss);
}

// --------------------------------------------------------------------------------------------- //

double getA(double theta) {
    return r/2 * cos(theta);
}

double getB(double theta) {
    return r/2 * sin(theta);
}

double getX(double theta1, double theta) {
    return getA(theta1) + r/2 * cos(theta);
}

double getY(double theta1, double theta) {
    return getB(theta1) + r/2 * sin(theta);
}

double pQuinaX(double x, double y) {
	double dis = dist(0,0,x,y);
	double h = sqrt((r/2)*(r/2) - (dis/2)*(dis/2));
	return x/2+h*y/dis;
}

double pQuinaY(double x, double y) {
	double dis = dist(0,0,x,y);
	double h = sqrt((r/2)*(r/2) - (dis/2)*(dis/2));
	return y/2-h*x/dis;
}

double anguloVetorEixoX(double x, double y) {

	if (y >= 0) return acos(x/(sqrt(x*x + y*y)));
	
	if (y < 0) return -acos(x/(sqrt(x*x + y*y)));

    return 0;
}

vector<point> removeDuplicates(vector<point> arr) {
    vector<point> unique;
    bool un = true;
    for(auto i:arr) {
        for(auto j:unique) {
            if (i.x == j.x && i.y == j.y) {
                un = false;
                break;
            }
        }

        if (un) unique.push_back(i);
        un = true;
    }

    return unique;
}

vector<point> normalizeCurve(vector<point> curva) {

    curva = removeDuplicates(curva);

    vector<point> newCurva;
    int s = curva.size();

    for (int i = 0; i < s; i++) {

        point pa = curva[i];
        point pp = curva[(i+1) % s];
        int qnt = floor(dist(pa.x,pa.y,pp.x,pp.y) / minDist);
		double dx = (pp.x - pa.x) / qnt;
		double dy = (pp.y - pa.y) / qnt;

        printf("criando pontos entre: (%.2f,%.2f) e (%.2f,%.2f)\n",pa.x,pa.y,pp.x,pp.y);

        for (int j = 0; j <= qnt; j++) newCurva.push_back({pa.x + dx * j, pa.y + dy * j});
    }

	curva = removeDuplicates(newCurva);
    cout << "----------------------------------------" << endl;
    cout << "Done. A curva esta com " << curva.size() << " pontos" << endl;
    cout << "----------------------------------------" << endl;
    return curva;
}

double calculaPonto(vector<point> curva, int idx, int s) {
    double u1,v1,u2,v2;
	point p0,p1,p2;

    p0 = curva[(idx-1 + s) % s];
    p1 = curva[(idx + s) % s];
    p2 = curva[(idx + 1 + s) % s];

    u1 = p2.x - p1.x;
    v1 = p2.y - p1.y;
    
    u2 = p1.x - p0.x;
    v2 = p1.y - p0.y;
	
	double cosTheta1 = -(c*u1+d*v1)/sqrt(u1*u1 + v1*v1);
    double cosTheta2 = -(c*u2+d*v2)/sqrt(u2*u2 + v2*v2);
	return (cosTheta1 + cosTheta2)*dist(p1.x,p1.y,p2.x,p2.y)*(r/2)/(2*sqrt(c*c + d*d));
}

void cGreen(vector<point> curva) {
	double res = 0;
    point p;
    int s = curva.size();

	for (int i = 0; i < s; i++) {
        p = curva[i];
		theta1 = anguloVetorEixoX(pQuinaX(p.x,p.y),pQuinaY(p.x,p.y));
		theta2 = anguloVetorEixoX(p.x-getA(theta1),p.y-getB(theta1));
		// calculando as coords
		a = getA(theta1);
        b = getB(theta1);
        x = getX(theta1,theta2);
        y = getY(theta1,theta2);
		c = -(y-b)/(r/2);
        d = (x-a)/(r/2);
		res += calculaPonto(curva, i, s);
	}

	resGreen = res;
}

void cGauss(vector<point> curva) {
	double res = 0;
    point p,q;
    int s = curva.size();

	for (int i = 0; i < curva.size(); i++) {
        p = curva[(i - 1 + s) % s];
        q = curva[(i + s) % s];
        res += p.x * q.y;
        res -= p.y * q.x;
	}

	resGauss = abs(res)/2;
}

double dist(double x1, double y1, double x2, double y2) {
    return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}